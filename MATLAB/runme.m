load JKS2e4

%Constants 
n_glen    = 3;
rele      = 1e-1;
eta_b     = 0.5;
eta_0     = 1.e+14/2.;
niter     = 5e6;
nout_iter = 1000;
epsi      = 3.171e-7;
damp = 0.96;
relaxation = 0.7;

%Initial guesses (except vx and vy that we already loaded)
etan   = 1.e+14*ones(nbe,1);
dVxdt  = zeros(nbv,1);
dVydt  = zeros(nbv,1);

%Manage derivatives once for all
%[areas alpha beta gamma]=NodalCoeffs(index,x,y,nbe);
weights=Weights(index,areas,nbe,nbv);

%Mesh size
[resolx resoly] = MeshSize(index,x,y,areas,weights,nbe,nbv);

%Physical properties once for all
[dsdx dsdy] = derive_xy_elem(surface,index,alpha,beta,nbe);
Helem       = mean(H(index),2);
rheology_B  = mean(rheology_B_temp(index),2);

%Update viscosity
[dvxdx dvxdy]=derive_xy_elem(vx,index,alpha,beta,nbe);
[dvydx dvydy]=derive_xy_elem(vy,index,alpha,beta,nbe);
for i=1:nbe
	eps_xx = dvxdx(i);
	eps_yy = dvydy(i);
	eps_xy = .5*(dvxdy(i)+dvydx(i));
	EII2 = eps_xx^2 + eps_yy^2 + eps_xy^2 + eps_xx*eps_yy;
	eta_it = 1.e+14/2.;
	if(EII2>0.) eta_it = rheology_B(i)/(2*EII2^((n_glen-1.)/(2*n_glen))); end
	etan(i) = min(eta_it,eta_0*1e5);
	if(isnan(etan(i))) error('Found NaN in etan(i)'); end
end

%Linear integration points order 3
wgt3=[0.555555555555556, 0.888888888888889, 0.555555555555556];
xg3 =[-0.774596669241483, 0.000000000000000, 0.774596669241483];

%Compute RHS amd ML once for all
ML            = zeros(nbv,1);
Fvx           = zeros(nbv,1);
Fvy           = zeros(nbv,1);
groundedratio = zeros(nbe,1);
isice         = false(nbe,1);
for n=1:nbe

	%Lumped mass matrix
	for i=1:3
		for j=1:3
			% \int_E phi_i * phi_i dE = A/6 and % \int_E phi_i * phi_j dE = A/12
			if(i==j)
				ML(index(n,j)) = ML(index(n,j))+areas(n)/6;
			else
				ML(index(n,j)) = ML(index(n,j))+areas(n)/12;
			end
		end
	end

	%Is there ice at all in the current element?
	level = ice_levelset(index(n,:));
	if level(1)<0 || level(2)<0 || level(3)<0
		isice(n)= true;
	else
		%We can skip this element altogether
		isice(n) = false;
		for i=1:3
			vx(index(n,i)) = 0.;
			vy(index(n,i)) = 0.;
		end
		continue;
	end

	%RHS (Driving stress term)
	for i=1:3
		for j=1:3
			if i==j
				Fvx(index(n,i)) = Fvx(index(n,i)) -rho*g*H(index(n,j))*dsdx(n)*areas(n)/6.;
				Fvy(index(n,i)) = Fvy(index(n,i)) -rho*g*H(index(n,j))*dsdy(n)*areas(n)/6.;
			else
				Fvx(index(n,i)) = Fvx(index(n,i)) -rho*g*H(index(n,j))*dsdx(n)*areas(n)/12.;
				Fvy(index(n,i)) = Fvy(index(n,i)) -rho*g*H(index(n,j))*dsdy(n)*areas(n)/12.;
			end
		end
	end

   %RHS (Water pressure at the ice front)
	level = ice_levelset(index(n,:));
	count = 0;
	for i=1:3
		if(level(i)<0.) count = count+1; end
	end
	if(count==1)
		%This element has an ice front, get indices of the 2 vertices
		seg1 = [index(n,1) index(n,2)];
		seg2 = [index(n,2) index(n,3)];
		seg3 = [index(n,3) index(n,1)];
		if     ice_levelset(seg1(1))>=0 && ice_levelset(seg1(2))>=0
			pairids = seg1;
		elseif ice_levelset(seg2(1))>=0 && ice_levelset(seg2(2))>=0
			pairids = seg2;
		elseif ice_levelset(seg3(1))>=0 && ice_levelset(seg3(2))>=0
			pairids = seg3;
		else
			error('not supported');
		end

		%Get normal
		len = sqrt((x(pairids(2))-x(pairids(1)))^2 + (y(pairids(2))-y(pairids(1)))^2);
		nx  = +(y(pairids(2))-y(pairids(1)))/len;
		ny  = -(x(pairids(2))-x(pairids(1)))/len;

		%add RHS
		for gg=1:3
			phi1 = (1.-xg3(gg))/2.;
			phi2 = (1.+xg3(gg))/2.;
			bg = base(pairids(1))*phi1 + base(pairids(2))*phi2;
			Hg = H(pairids(1))*phi1 + H(pairids(2))*phi2;
			bg = min(bg,0);
			Fvx(pairids(1)) = Fvx(pairids(1)) +wgt3(gg)/2*1/2*(-rho_w*g*bg^2+rho*g*Hg^2)*nx*len*phi1;
			Fvx(pairids(2)) = Fvx(pairids(2)) +wgt3(gg)/2*1/2*(-rho_w*g*bg^2+rho*g*Hg^2)*nx*len*phi2;
			Fvy(pairids(1)) = Fvy(pairids(1)) +wgt3(gg)/2*1/2*(-rho_w*g*bg^2+rho*g*Hg^2)*ny*len*phi1;
			Fvy(pairids(2)) = Fvy(pairids(2)) +wgt3(gg)/2*1/2*(-rho_w*g*bg^2+rho*g*Hg^2)*ny*len*phi2;
		end
	end

	%One more thing in this element loop: prepare groundedarea needed later for the calculation of basal friction
	level = ocean_levelset(index(n,:));
	if level(1)>=0 && level(2)>=0 && level(3)>=0
		%Completely grounded
		groundedratio(n)=1.;
	elseif level(1)<=0 && level(2)<=0 && level(3)<=0
		%Completely floating
		groundedratio(n)=0.;
	else
		%Partially floating,
		if(level(1)*level(2)>0) %Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2
			s1=level(3)/(level(3)-level(2));
			s2=level(3)/(level(3)-level(1));
		elseif(level(2)*level(3)>0) %Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2
			s1=level(1)/(level(1)-level(2));
			s2=level(1)/(level(1)-level(3));
		elseif(level(1)*level(3)>0) %Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2
			s1=level(2)/(level(2)-level(1));
			s2=level(2)/(level(2)-level(3));
		else
			error('not supposed to be here...');
		end

		if(level(1)*level(2)*level(3)>0)
			groundedratio(n)= s1*s2;
		else
			groundedratio(n)= (1-s1*s2);
		end
	end
end

%Finally add calculation of friction coefficient:
alpha2 = zeros(nbv,1);
for i=1:nbv

   %Compute effective pressure
   p_ice   = g*rho*H(i);
   p_water = -rho_w*g*base(i);
   Neff    = p_ice - p_water;
   if(Neff<0.) Neff=0.; end

   %Compute alpha2
   alpha2(i) = friction(i)^2*Neff;
end


%%Modifications%%
head = zeros(nbv,1);
next = zeros(nbe*3, 1);


for i = 1:nbv
     head(i) = -1;
end

%Now construct the chain
for k=1:nbe
    for j=1:3
        i = index(k,j);  
        p = 3*(k-1)+j;
        next(p) = head(i);
        head(i) = p;       
    %    fprintf('i=%d, next=%d, head=%d \n',i, next(p), head(i));
    end
end

%Now we can construct the connectivity matrix
MAXCONNECT = 8;
connectivity = zeros(nbv*MAXCONNECT, 1);
columns = zeros(nbv*MAXCONNECT, 1);
for i=1:nbv

   %Go over all of the elements connected to node I%
   count = 1;
   p = head(i);
  while (p~=-1)
   
      k = fix(p/3) + 1;     %”row" in index
      j = rem(p, 3);   %"column" in index
      
      if (j==0)
          j=3;
          k= k -1;
      end
      %sanity check
      assert(index(k,j)==i);
      
       	%fprintf('k=%d, j=%d, i=%d \n',k, j, i);
      %enter element in connectivity matrix
      connectivity((i-1)* MAXCONNECT + count) = k;
       columns((i-1)* MAXCONNECT + count) = j;
      count = count + 1;
     p = next(p); 
   end
  
end

%Main loop, allocate a few vectors needed for the computation

tic
for iter = 1:niter % Pseudo-Transient cycles

	%Strain rates
	[dvxdx dvxdy]=derive_xy_elem(vx,index,alpha,beta,nbe);
	[dvydx dvydy]=derive_xy_elem(vy,index,alpha,beta,nbe);

	%KV term in equation 22
	KVx  = zeros(nbv,1);
	KVy  = zeros(nbv,1);
	kvx =  zeros(nbe,3);
    kvy =  zeros(nbe,3);

    for n=1:nbe

		%Skip if no ice
		if ~isice(n) continue; end

		%Viscous Deformation
		eta_e  = etan(n);
		eps_xx = dvxdx(n);
		eps_yy = dvydy(n);
		eps_xy = .5*(dvxdy(n)+dvydx(n));
		for i=1:3
			kvx(n,i) = 2*Helem(n)*eta_e*(2*eps_xx+eps_yy)*alpha(n,i)*areas(n) + 2*Helem(n)*eta_e*eps_xy*beta(n,i)*areas(n);
            kvy(n,i) = 2*Helem(n)*eta_e*eps_xy*alpha(n,i)*areas(n) + 2*Helem(n)*eta_e*(2*eps_yy+eps_xx)*beta(n,i)*areas(n);
        end
        
        %Add basal friction
		if groundedratio(n)>0.
         %n3 = n * 3;
         gr_a = groundedratio(n) * areas(n);

                for k=1:3
                    for i = 1:3
                         i_index = index(n,i);
                         gr_a_alpha2 = gr_a*alpha2(i_index);

                        for j=1:3 
                            j_index = index(n,j);
                            gr_a_alpha2_vx = gr_a_alpha2 * vx(j_index);
                            gr_a_alpha2_vy = gr_a_alpha2 * vy(j_index);
                            if (i == j && j == k) 
                                kvx(n, k) = kvx(n, k) + gr_a_alpha2_vx / 10.;
                                kvy(n, k) = kvy(n, k) + gr_a_alpha2_vy / 10.;
                            elseif ((i ~= j) && (j ~= k) && (k ~= i)) 
                                kvx(n, k) = kvx(n, k) + gr_a_alpha2_vx / 60.;
                                kvy(n, k) = kvy(n, k) + gr_a_alpha2_vy / 60.;
                            else 
                                kvx(n, k) = kvx(n, k) + gr_a_alpha2_vx / 30.;
                                kvy(n, k) = kvy(n, k) + gr_a_alpha2_vy / 30.;
                                end

                            end

                        end
                    end
               end          
        end
    

    for i = 1:nbv
        for j = 1:8
            if connectivity((i-1)*8+j) ~=0 
                   KVx(i) = KVx(i) + kvx((connectivity((i-1)*8+j)), (columns((i-1)*8+j)));
                   KVy(i) = KVy(i) + kvy((connectivity((i-1)*8+j)), (columns((i-1)*8+j)));
            end
        end
    end
   
    
	%Get current viscosity on nodes (Needed for time stepping)
    eta_nbv = zeros(nbv,1);
    Eta_nbv = zeros(nbe,1);
    
    for i = 1:nbe
            Eta_nbv(i) = etan(i)*areas(i);
    end
    
    
        for i = 1:nbv
        for j = 1:8
            if connectivity((i-1)*8+j) ~=0 
            eta_nbv(i) = eta_nbv(i) + Eta_nbv(connectivity((i-1)*8+j));  
            end
        end
    end
    
    
    for i=1:nbv
        eta_nbv(i) = eta_nbv(i)/weights(i);
    end

	%Velocity rate update in the x and y, refer to equation 19 in Rass paper
	normX = 0;
	normY = 0;

	%1. Get time derivative based on residual (dV/dt)
	ResVx =  1./(rho*ML.*max(80,H)).*(-KVx + Fvx); %rate of velocity in the x, equation 23
	ResVy =  1./(rho*ML.*max(80,H)).*(-KVy + Fvy); %rate of velocity in the y, equation 24
	

    dVxdt = dVxdt*damp + ResVx;
	dVydt = dVydt*damp + ResVy;
    
	if(any(isnan(dVxdt))) error('Found NaN in dVxdt[i]'); end
	if(any(isnan(dVydt))) error('Found NaN in dVydt[i]'); end

	%2. Explicit CFL time step for viscous flow, x and y directions
    dtVx = rho*resolx.^2./(4*eta_nbv*(1.+eta_b)*4.1);
    dtVy = rho*resoly.^2./(4*eta_nbv*(1.+eta_b)*4.1);
    
	%3. velocity update, vx(new) = vx(old) + change in vx, Similarly for vy
	vx = vx + relaxation*dVxdt.*dtVx;
	vy = vy + relaxation*dVydt.*dtVy;

	%Apply Dirichlet boundary condition, Residual should also be 0 (for convergence)
	pos = find(~isnan(spcvx));
	vx(pos) = spcvx(pos);
	dVxdt(pos) = 0.;
	pos = find(~isnan(spcvy));
	vy(pos) = spcvy(pos);
	dVydt(pos) = 0.;


	%LAST: Update viscosity
	for i=1:nbe
		if ~isice(i) continue; end
		eps_xx = dvxdx(i);
		eps_yy = dvydy(i);
		eps_xy = .5*(dvxdy(i)+dvydx(i));
		EII2 = eps_xx^2 + eps_yy^2 + eps_xy^2 + eps_xx*eps_yy;
		eta_it = 1.e+14/2.;
		if(EII2>0.) eta_it = rheology_B(i)/(2*EII2^((n_glen-1.)/(2*n_glen))); end

		etan(i) = min(exp(rele*log(eta_it) + (1-rele)*log(etan(i))),eta_0*1e5);
		if(isnan(etan(i))) error('Found NaN in etan(i)'); end
    end
      %  plotmodel(md,'data',sqrt(vx.^2+vy.^2)*yts,'data', dtVx);
        iterror = max(max(abs(dVxdt)),max(abs(dVydt)));
   
	if((iterror < epsi) && (iter > 2)) break; end
	if (mod(iter,nout_iter)==1)
		fprintf('iter=%d, err=%1.3e \n',iter,iterror)
		clf
		%plotmodel(md,'data',sqrt(vx.^2+vy.^2)*yts); %Multiplication by yts converts m/s to m/yr
		drawnow
    end   
end

clf
%plotmodel(md,'data',sqrt(vx.^2+vy.^2)*yts);
fprintf('iter=%d, err=%1.3e --> converged\n',iter,iterror)
toc

 

function weights=Weights(index,areas,nbe,nbv)% {{{
	weights = sparse(index(:),ones(3*nbe,1),repmat(areas,3,1),nbv,1);
end % }}}
function [dfdx_e dfdy_e] = derive_xy_elem(f,index,alpha,beta,nbe)% {{{

	dfdx_e = sum(f(index).*alpha,2);
	dfdy_e = sum(f(index).*beta,2);

end % }}}
function f_v = elem2node(f_e,index,areas,weights,nbe,nbv)% {{{

	f_v = sparse(index(:),ones(3*nbe,1),repmat(areas.*f_e,3,1),nbv,1)./weights;

end % }}}
function [resolx resoly] = MeshSize(index,x,y,areas,weights,nbe,nbv)% {{{

	%Get element size
	dx_elem=max(x(index),[],2)-min(x(index),[],2);
	dy_elem=max(y(index),[],2)-min(y(index),[],2);

	%Average over each node
	resolx = elem2node(dx_elem,index,areas,weights,nbe,nbv);
	resoly = elem2node(dy_elem,index,areas,weights,nbe,nbv);
end % }}}
