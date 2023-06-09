#include <iostream>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>
using namespace std;

/*I/O stuff {{{*/
FILE* SetFilePointerToData(FILE* fid,int* pcode,int* pvector_type,const char* data_name){/*{{{*/

    int found  = 0;
    const char* mddot = "md.";
    char* record_name = NULL;
    int   record_name_size;
    long long record_length;
    int record_code;       //1 to 7 number
    int vector_type   = 0; //nodal or elementary

    if(strncmp(data_name,mddot,3)!=0){
        std::cerr <<"Cannot fetch \""<<data_name<<"\" does not start with \""<<mddot<<"\"";
    }

    /*First set FILE* position to the beginning of the file: */
    fseek(fid,0,SEEK_SET);

    /*Now march through file looking for the correct data identifier: */
    for(;;){
        /*Read size of first string name: */
        if(fread(&record_name_size,sizeof(int),1,fid)==0){
            /*we have reached the end of the file. break: */
            delete record_name;
            break;
        }
        if(record_name_size<3 || record_name_size>80){
            std::cerr<<"error while looking in binary file. Found a string of size "<<record_name_size;
        }

        /*Allocate string of correct size: */
        record_name=new char[record_name_size+1];
        record_name[record_name_size]='\0';

        /*Read record_name: */
        if(fread(record_name,record_name_size*sizeof(char),1,fid)==0){
            /*we have reached the end of the file. break: */
            found=0;
            delete [] record_name;
            break;
        }
        if(strncmp(record_name,mddot,3)!=0){
            std::cerr<<"error while reading binary file: record does not start with \"md.\": "<<record_name;
        }

        /*Is this the record sought for? : */
        if(strcmp(record_name,data_name)==0){
            /*Ok, we have found the correct string. Pass the record length, and read data type code: */
            fseek(fid,sizeof(long long),SEEK_CUR);
            if(fread(&record_code,sizeof(int),1,fid)!=1) std::cerr<<"Could not read record_code";
            /*if record_code points to a vector, get its type (nodal or elementary): */
            if((5<=record_code && record_code<=7) || record_code==10){
                if(fread(&vector_type,sizeof(int),1,fid)!=1) std::cerr<<"Could not read vector_type";
            }
            found=1;
            delete [] record_name;
            break;
        }
        else{
            /*This is not the correct string, read the record length, and use it to skip this record: */
            if(fread(&record_length,sizeof(long long),1,fid)!=1) std::cerr<<"Could not read record_length";
            /*skip: */
            fseek(fid,record_length,SEEK_CUR);
            delete [] record_name;
        }
    }
    if(!found) std::cerr<<"could not find data with name \"" << data_name << "\" in binary file";

    /*Assign output pointers:*/
    *pcode=record_code;
    if(pvector_type) *pvector_type=vector_type;

    return fid;
}
/*}}}*/
void  FetchData(FILE* fid,int* pinteger,const char* data_name){/*{{{*/

    /*output: */
    int integer;
    int code;

    /*Set file pointer to beginning of the data: */
    fid=SetFilePointerToData(fid,&code,NULL,data_name);

    if(code!=2)std::cerr <<"expecting an integer for \"" << data_name<<"\"";

    /*We have to read a integer from disk. First read the dimensions of the integer, then the integer: */
    if(fread(&integer,sizeof(int),1,fid)!=1) std::cerr<<"could not read integer ";

    /*Assign output pointers: */
    *pinteger=integer;
}/*}}}*/
void  FetchData(FILE* fid,int** pmatrix,int* pM,int* pN,const char* data_name){/*{{{*/

    /*output: */
    int M,N;
    double* matrix=NULL;
    int* integer_matrix=NULL;
    int code=0;

    /*Set file pointer to beginning of the data: */
    fid=SetFilePointerToData(fid,&code,NULL,data_name);
    if(code!=5 && code!=6 && code!=7)std::cerr<<"expecting a IssmDouble, integer or boolean matrix for \""<<data_name<<"\""<<" (Code is "<<code<<")";

    /*Now fetch: */

    /*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */
    /*numberofelements: */
    if(fread(&M,sizeof(int),1,fid)!=1) std::cerr<<"could not read number of rows for matrix ";
    if(fread(&N,sizeof(int),1,fid)!=1) std::cerr<<"could not read number of columns for matrix ";

    /*Now allocate matrix: */
    if(M*N){
        matrix=new double[M*N];

        /*Read matrix on node 0, then broadcast: */
        if(fread(matrix,M*N*sizeof(double),1,fid)!=1) std::cerr<<"could not read matrix ";
    }

    /*Now cast to integer: */
    if(M*N){
        integer_matrix=new int[M*N];
        for (int i=0;i<M;i++){
            for (int j=0;j<N;j++){
                integer_matrix[i*N+j]=(int)matrix[i*N+j];
            }
        }
    }
    else{
        integer_matrix=NULL;
    }
    /*Free ressources:*/
    delete [] matrix;

    /*Assign output pointers: */
    *pmatrix=integer_matrix;
    if(pM)*pM=M;
    if(pN)*pN=N;
}/*}}}*/
void  FetchData(FILE* fid,double* pdouble,const char* data_name){/*{{{*/

    /*output: */
    double value;
    int code;

    /*Set file pointer to beginning of the data: */
    fid=SetFilePointerToData(fid,&code,NULL,data_name);

    if(code!=3)std::cerr <<"expecting a double for \"" << data_name<<"\"";

    /*We have to read a integer from disk. First read the dimensions of the integer, then the integer: */
    if(fread(&value,sizeof(double),1,fid)!=1) std::cerr<<"could not read scalar";

    /*Assign output pointers: */
    *pdouble=value;
}/*}}}*/
void  FetchData(FILE* fid,double** pmatrix,int* pM,int* pN,const char* data_name){/*{{{*/

    /*output: */
    int M,N;
    double* matrix=NULL;
    int* integer_matrix=NULL;
    int code=0;

    /*Set file pointer to beginning of the data: */
    fid=SetFilePointerToData(fid,&code,NULL,data_name);
    if(code!=5 && code!=6 && code!=7)std::cerr<<"expecting a IssmDouble, integer or boolean matrix for \""<<data_name<<"\""<<" (Code is "<<code<<")";

    /*Now fetch: */

    /*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */
    /*numberofelements: */
    if(fread(&M,sizeof(int),1,fid)!=1) std::cerr<<"could not read number of rows for matrix ";
    if(fread(&N,sizeof(int),1,fid)!=1) std::cerr<<"could not read number of columns for matrix ";

    /*Now allocate matrix: */
    if(M*N){
        matrix=new double[M*N];

        /*Read matrix on node 0, then broadcast: */
        if(fread(matrix,M*N*sizeof(double),1,fid)!=1) std::cerr<<"could not read matrix ";
    }

    /*Assign output pointers: */
    *pmatrix=matrix;
    if(pM)*pM=M;
    if(pN)*pN=N;
}/*}}}*/
void  WriteData(FILE* fid,double* matrix,int M,int N,const char* data_name){/*{{{*/

    /*First write enum: */
    int length=(strlen(data_name)+1)*sizeof(char);
    fwrite(&length,sizeof(int),1,fid);
    fwrite(data_name,length,1,fid);

    /*Now write time and step: */
    double time = 0.;
    int    step = 1;
    fwrite(&time,sizeof(double),1,fid);
    fwrite(&step,sizeof(int),1,fid);

    /*writing a IssmDouble array, type is 3:*/
    int type=3;
    fwrite(&type,sizeof(int),1,fid);
    fwrite(&M,sizeof(int),1,fid);
    fwrite(&N,sizeof(int),1,fid);
    fwrite(matrix,M*N*sizeof(double),1,fid);
}/*}}}*/
void  WriteData(FILE* fid,const char* string,const char* data_name){/*{{{*/

    /*First write enum: */
    int length=(strlen(data_name)+1)*sizeof(char);
    fwrite(&length,sizeof(int),1,fid);
    fwrite(data_name,length,1,fid);

    /*Now write time and step: */
    double time = 0.;
    int    step = 1;
    fwrite(&time,sizeof(double),1,fid);
    fwrite(&step,sizeof(int),1,fid);

    /*writing a string, type is 2: */
    int type=2;
    fwrite(&type,sizeof(int),1,fid);

    length=(strlen(string)+1)*sizeof(char);
    fwrite(&length,sizeof(int),1,fid);
    fwrite(string,length,1,fid);
}/*}}}*/
/*}}}*/
void NodalCoeffs(double** pareas,double** palpha,double** pbeta,int* index,double* x,double* y,int nbe){/*{{{*/

    /*Allocate output vectors*/
    double* areas = new double[nbe];
    double* alpha = new double[nbe*3];
    double* beta  = new double[nbe*3];
    for(int i = 0; i < nbe; i++) {
        int n1 = index[i*3+0]-1;
        int n2 = index[i*3+1]-1;
        int n3 = index[i*3+2]-1;

        double x1 = x[n1];
        double x2 = x[n2];
        double x3 = x[n3];
        double y1 = y[n1];
        double y2 = y[n2];
        double y3 = y[n3];

        double invdet = 1./(x1 * (y2 - y3) - x2 * (y1 - y3) + x3 * (y1 - y2));

        alpha[i*3+0] = invdet * (y2 - y3);
        alpha[i*3+1] = invdet * (y3 - y1);
        alpha[i*3+2] = invdet * (y1 - y2);

        beta[i*3+0] = invdet * (x3 - x2);
        beta[i*3+1] = invdet * (x1 - x3);
        beta[i*3+2] = invdet * (x2 - x1);

        areas[i]= 0.5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1));
    }

    /*Assign output pointers*/
    *pareas = areas;
    *palpha = alpha;
    *pbeta = beta;
}/*}}}*/
void Weights(double** pweights,int* index,double* areas,int nbe,int nbv){/*{{{*/

    /*Allocate output and initialize as 0*/
    double* weights = new double[nbv];
    for(int i = 0; i < nbv; i++) weights[i]=0.;

    /*Loop over elements*/
    for(int i = 0; i < nbe; i++){
        for(int j = 0; j < 3; j++){
            weights[index[i*3+j]-1] += areas[i];
        }
    }

    /*Assign output pointer*/
    *pweights = weights;
}/*}}}*/
void derive_xy_elem(double* dfdx_e,double* dfdy_e,double* f,int* index,double* alpha,double* beta,int nbe){/*{{{*/

    /*WARNING!! Assume that dfdx_e and dfdy_e have been properly allocated*/
    for(int i=0;i<nbe;i++){
        int n1 = index[i*3+0]-1;
        int n2 = index[i*3+1]-1;
        int n3 = index[i*3+2]-1;
        dfdx_e[i] = f[n1]*alpha[i*3+0] + f[n2]*alpha[i*3+1] + f[n3]*alpha[i*3+2];
        dfdy_e[i] = f[n1]*beta[ i*3+0] + f[n2]*beta[ i*3+1] + f[n3]*beta[ i*3+2];
    }
}/*}}}*/
void elem2node(double* f_v,double* f_e,int* index,double* areas,double* weights,int nbe,int nbv){/*{{{*/

    /*WARNING!! Assume that f_v has been properly allocated*/

    /*Reinitialize output*/
    for(int i=0;i<nbv;i++) f_v[i] = 0.;

    /*Add contributions from all elements connected to vertex i*/
    for(int i=0;i<nbe;i++){
        int n1 = index[i*3+0]-1;
        int n2 = index[i*3+1]-1;
        int n3 = index[i*3+2]-1;
        f_v[n1] += f_e[i]*areas[i];
        f_v[n2] += f_e[i]*areas[i];
        f_v[n3] += f_e[i]*areas[i];
    }

    /*Divide by sum of areas*/
    for(int i=0;i<nbv;i++) f_v[i] = f_v[i]/weights[i];

}/*}}}*/
void MeshSize(double* resolx,double* resoly,int* index,double* x,double* y,double* areas,double* weights,int nbe,int nbv){/*{{{*/

    /*Get element size*/
    double  xmin,xmax,ymin,ymax;
    double* dx_elem = new double[nbe];
    double* dy_elem = new double[nbe];
    for(int i=0;i<nbe;i++){
        int n1 = index[i*3+0]-1;
        int n2 = index[i*3+1]-1;
        int n3 = index[i*3+2]-1;
        xmin = min(min(x[n1],x[n2]),x[n3]);
        xmax = max(max(x[n1],x[n2]),x[n3]);
        ymin = min(min(y[n1],y[n2]),y[n3]);
        ymax = max(max(y[n1],y[n2]),y[n3]);
        dx_elem[i] = xmax - xmin;
        dy_elem[i] = ymax - ymin;
    }

    /*Average over each node*/
    elem2node(resolx,dx_elem,index,areas,weights,nbe,nbv);
    elem2node(resoly,dy_elem,index,areas,weights,nbe,nbv);

    /*Cleanup and return*/
    delete [] dx_elem;
    delete [] dy_elem;
}/*}}}*/

// timer
#include "sys/time.h"
double timer_start = 0;
double cpu_sec(){ struct timeval tp; gettimeofday(&tp,NULL); return tp.tv_sec+1e-6*tp.tv_usec; }
void   tic(){ timer_start = cpu_sec(); }
double toc(){ return cpu_sec()-timer_start; }
void   tim(const char *what, double n){ double s=toc(); printf("%s: %8.3f seconds",what,s);if(n>0)printf(", %8.3f GB/s", n/s); printf("\n"); }

void meanForHAndRheologyB(double* outputArray, double* inputArray, int* index, int arraySize)
{
    for (int i=0; i<arraySize; i++)
    {
        double value = 0;
        for (int j=0; j<arraySize; j++)
        {
            value += inputArray[index[i*3+j]];
        }
        double mean = (double)(value/3.0);
        outputArray[i] = mean;
    }
}
int main() {/*{{{*/

    /*Open input binary file*/
    const char *inputfile = "/home/caelinux2/CLionProjects/bin_files/JKS2e4.bin";
    const char *outputfile = "/home/caelinux2/CLionProjects/V2/output.outbin";
    FILE *fid = fopen(inputfile, "rb");
    if (fid == NULL) std::cerr << "could not open file " << inputfile << " for binary reading or writing";


    /*Get All we need from binary file*/
    int nbe, nbv, M, N;
    double g, rho, rho_w, yts;
    int *index = NULL;
    double *spcvx = NULL;
    double *spcvy = NULL;
    double *x = NULL;
    double *y = NULL;
    double *H = NULL;
    double *surface = NULL;
    double *base = NULL;
    double *ice_levelset = NULL;
    double *ocean_levelset = NULL;
    double *rheology_B_temp = NULL;
    double *vx = NULL;
    double *vy = NULL;
    double *friction = NULL;
    FetchData(fid, &nbe, "md.mesh.numberofelements");
    FetchData(fid, &nbv, "md.mesh.numberofvertices");
    FetchData(fid, &g, "md.constants.g");
    FetchData(fid, &rho, "md.materials.rho_ice");
    FetchData(fid, &rho_w, "md.materials.rho_water");
    FetchData(fid, &yts, "md.constants.yts");
    FetchData(fid, &index, &M, &N, "md.mesh.elements");
    FetchData(fid, &spcvx, &M, &N, "md.stressbalance.spcvx");
    FetchData(fid, &spcvy, &M, &N, "md.stressbalance.spcvy");
    FetchData(fid, &x, &M, &N, "md.mesh.x");
    FetchData(fid, &y, &M, &N, "md.mesh.y");
    FetchData(fid, &H, &M, &N, "md.geometry.thickness");
    FetchData(fid, &surface, &M, &N, "md.geometry.surface");
    FetchData(fid, &base, &M, &N, "md.geometry.base");
    FetchData(fid, &ice_levelset, &M, &N, "md.mask.ice_levelset");
    FetchData(fid, &ocean_levelset, &M, &N, "md.mask.ocean_levelset");
    FetchData(fid, &rheology_B_temp, &M, &N, "md.materials.rheology_B");
    FetchData(fid, &vx, &M, &N, "md.initialization.vx");
    FetchData(fid, &vy, &M, &N, "md.initialization.vy");
    FetchData(fid, &friction, &M, &N, "md.friction.coefficient");

    /*Close input file*/
    if (fclose(fid) != 0) std::cerr << "could not close file " << inputfile;


    /*Constants*/
    double n_glen = 3.;
    double damp = 0.96;  //0.96
    double rele = 1e-1;
    double eta_b = 0.5;
    double eta_0 = 1.e+14 / 2.;
    int niter = 5e6;  //5e6
    int nout_iter = 100;
    double epsi = 3.171e-7;
    double relaxation = 0.7;

    /*Initial guesses (except vx and vy that we already loaded)*/
    double *etan = new double[nbe];
    for (int i = 0; i < nbe; i++) etan[i] = 1.e+14;
    double *dVxdt = new double[nbv];
    for (int i = 0; i < nbv; i++) dVxdt[i] = 0.;
    double *dVydt = new double[nbv];
    for (int i = 0; i < nbv; i++) dVydt[i] = 0.;

    /*Manage derivatives once for all*/
    double *alpha = NULL;
    double *beta = NULL;
    double *areas = NULL;
    double *weights = NULL;
    NodalCoeffs(&areas, &alpha, &beta, index, x, y, nbe);
    Weights(&weights, index, areas, nbe, nbv);


    /*MeshSize*/
    double *resolx = new double[nbv];
    double *resoly = new double[nbv];
    MeshSize(resolx, resoly, index, x, y, areas, weights, nbe, nbv);

    //Marco's code
    double* dsdx = new double[nbe];
    double* dsdy = new double[nbe];
    derive_xy_elem(dsdx, dsdy, surface, index, alpha, beta, nbe);
    double* Helem = new double[nbv];
    meanForHAndRheologyB(Helem, H, index, nbv);// = mean(H(index), 2);
    double *rheology_B = new double[nbv];
    meanForHAndRheologyB(rheology_B, rheology_B_temp, index, nbv);// = mean(rheology_B_temp(index), 2) //no clue how this works
    //update viscosity
    double* dvxdx = new double[nbe];
    double* dvydy = new double[nbe];
    double* dvxdy = new double[nbe];
    double* dvydx = new double[nbe];
    derive_xy_elem(dvxdx, dvxdy, vx, index, alpha, beta, nbe);//[dvxdx dvxdy]=derive_xy_elem(vx,index,alpha,beta,nbe); MATLAB CODE
    derive_xy_elem(dvydx, dvydy, vy, index, alpha, beta, nbe);//[dvydx dvydy]=derive_xy_elem(vy,index,alpha,beta,nbe);

    for (int i=0; i<nbe; i++)
    {
        double eps_xx = dvxdx[i];
        double eps_yy = dvydy[i];
        double eps_xy = .5*(dvxdy[i] + dvydx[i]);
        double EII2 = eps_xx*eps_xx + eps_yy*eps_yy + eps_xy*eps_xy + eps_xx*eps_yy;
        double eta_it = 50000000000000; //should be 1.0e+14/2
        if (EII2 > 0)
        {
            eta_it = rheology_B[i]/(2*pow(EII2, ((n_glen-1))/(2*n_glen))); //rheology_B doesnt exist?
        }
        etan[i] = std::min(eta_it, eta_0*100000); //10000 is supposed to be 1e5
        if (isnan(etan[i])) std::cout << "found NaN in etan[i]\n"; //if(isnan(etan(i))) error('Found NaN in etan(i)'); end
    }

    //Linear integration points order 3
    double wgt3[3] = {0.555555555555556, 0.888888888888889, 0.555555555555556};
    double xg3[3] = {-0.774596669241483, 0.000000000000000, 0.774596669241483};


    //Compute RHS and ML once for all
    double* ML = new double[nbv];
    double* Fvx = new double[nbv];
    double* Fvy = new double[nbv];
    double* groundedratio = new double[nbe];
    double* isice = new double[nbe];

    for (int n=0; n<nbe; n++)
    {
        //Lumped mass matrix
        for (int i=0; i<3; i++)
        {
            for (int j=1; j<3; j++)
            {
                if (i==j) ML[index[n*3+j]] = ML[index[n*3+j]] + areas[n]/6;//ML(index(n,j)) = ML(index(n,j))+areas(n)/6;
                else ML[index[n*3+j]] = ML[index[n*3+j]] + areas[n]/12;//ML(index(n,j)) = ML(index(n,j))+areas(n)/12
            }
        }
        //is there ice at all in the current element
        double* level = new double[3];
        for (int i=0; i<(3); i++) {level[i] = ice_levelset[index[n*3+i]];} // level = ice_levelset(index(n, :));
        if (level[0]<0 or level[1]<0 or level[2] <0) isice[n] = true;
        else
        {
            //we can skip this element altogether
            isice[n] = false;
            for (int i=0; i<3; i++)
            {
                vx[index[n*3+i]] = 0.;
                vy[index[n*3+i]] = 0.;
            }
            continue;
        }

        //RHS (driving stress term)
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {
                if (i==j)
                {
                    Fvx[index[n*3+i]] = Fvx[index[n*3+i]] - rho*g*H[index[n*3+j]] * dsdx[n] * areas[n]/6.0;
                    Fvy[index[n*3+i]] = Fvy[index[n*3+i]] - rho*g*H[index[n*3+j]] * dsdy[n] * areas[n]/6.0;
                }
                else
                {
                    Fvx[index[n*3+i]] = Fvx[index[n*3+i]] - rho*g*H[index[n*3+j]] * dsdx[n] * areas[n]/12.0;
                    Fvy[index[n*3+i]] = Fvy[index[n*3+i]] - rho*g*H[index[n*3+j]] * dsdy[n] * areas[n]/12.0;
                }
            }
        }

        //RHS (water pressure at the ice front)
        for (int i=0; i<(3); i++) {level[i] = ice_levelset[index[n*3+i]];} // level = ice_levelset(index(n, :));
        int count = 0;
        for (int i=0; i<3; i++)
        {
            if (level[i]<0.0) count = count+1;
        }

        //define pairids here
        int* pairids;
        if (count==1)
        {
            //this element has an ice front, get indices of the 2 vertices
            int seg1[2] = {index[n*3+0], index[n*3+1]};//[index(n, 1) index(n,2)];
            int seg2[2] = {index[n*3+1], index[n*3+2]};//[index(n, 2) index(n,3)];
            int seg3[2] = {index[n*3+2], index[n*3+0]};//[index(n, 3) index(n,1)];
            if (ice_levelset[seg1[0]] >=0 and ice_levelset[seg1[1]]>=0) pairids = seg1;
            else if (ice_levelset[seg2[0]] >=0 and ice_levelset[seg2[1]]>=0) pairids = seg2;
            else if (ice_levelset[seg3[0]] >=0 and ice_levelset[seg3[1]]>=0) pairids = seg3;
            else std::cout << "not supported\n";

            //get normal
            double len = sqrt(pow((x[(int)pairids[1]] - x[(int)pairids[0]]), 2) + pow((y[(int)pairids[1]] - y[(int)pairids[0]]), 2));
            double nx = +(y[(int)pairids[1]]-y[(int)pairids[0]])/len;
            double ny = -(x[(int)pairids[1]]-x[(int)pairids[0]])/len;

            //add RHS
            for (int gg=0; gg<3; gg++)
            {
                double phi1 = (1.0-xg3[gg])/2.0;
                double phi2 = (1.0+xg3[gg])/2.0;
                double bg = base[(int)pairids[0]] * phi1 + base[(int)pairids[1]] * phi2;
                double Hg = H[(int)pairids[0]] * phi1 + H[(int)pairids[1]] * phi2;
                //bg; //min(bg, 0);
                Fvx[(int)pairids[0]] = Fvx[(int)pairids[0]] + wgt3[gg]/2.0*1/2.0*(pow(-rho_w*g*bg, 2)+pow(rho*g*Hg, 2))*nx*len*phi1;
                Fvx[(int)pairids[1]] = Fvx[(int)pairids[1]] + wgt3[gg]/2.0*1/2.0*(pow(-rho_w*g*bg, 2)+pow(rho*g*Hg, 2))*nx*len*phi2;
                Fvy[(int)pairids[0]] = Fvx[(int)pairids[0]] + wgt3[gg]/2.0*1/2.0*(pow(-rho_w*g*bg, 2)+pow(rho*g*Hg, 2))*nx*len*phi1;
                Fvy[(int)pairids[1]] = Fvx[(int)pairids[1]] + wgt3[gg]/2.0*1/2.0*(pow(-rho_w*g*bg, 2)+pow(rho*g*Hg, 2))*nx*len*phi2;
            }
        }
        //One more thing in this element loop: prepare groundedarea needed later for the calculations of basal friction
        for (int i=0; i<(3); i++) {level[i] = ocean_levelset[index[n*3+i]];} // level = ocean_levelset[index[n, :]];
        if (level[0] >=0 and level[1]>=0 and level[2] >=0)
        {
            //completely grounded
            groundedratio[n] = 1.0;
        }
        else if (level[0]<=0 and level[1]<=0 and level[2] <=0)
        {
            //Completely floating
            groundedratio[n] = 0.0;
        }
        else
        {
            //partially floating,
            double s1,s2;
            if (level[0]*level[1]>0) //Nodes 0 and 1 are similiar, so points must be found on segment 0-2 and 1-2
            {
                s1 = level[2]/(level[2]-level[1]);
                s2 = level[2]/(level[2]-level[0]);
            }
            else if (level[1]*level[2]>0) //Nodes 1 and 2 are similiar, so points must be found on segment 0-1 ad 0-2
            {
                s1 = level[0]/(level[0]-level[1]);
                s2 = level[0]/(level[0]-level[2]);
            }
            else if (level[0]*level[2]>0) //Nodes 0 and 2 are similiar, so points must be found on segment 1-0 and 1-2
            {
                s1 = level[1]/(level[1]-level[0]);
                s2 = level[1]/(level[1]-level[2]);
            }
            else std::cout << "Not Supposed to be here...\n";

            if (level[0]*level[1]*level[2]>0) groundedratio[n] = s1*s2;
            else groundedratio[n] = (1-s1*s2);
        }
    }

    /* Kevins Portion */
    /* Finally add claculations of friction coefficient: */
    double *alpha2 = new double[nbv];
    for (int i = 0; i<nbv; i++) {

        /*compute effective pressure*/
        double p_ice = g*rho*H[i];
        double p_water = -rho_w*g*base[i];
        double Neff = p_ice - p_water;
        if (Neff < 0) Neff = 0;

        /*compute alpha2*/
        alpha2[i] = pow(friction[i], 2)*Neff;
    }

    /* Modifications */
    int *head = new int[nbv];
    int *next = new int[(nbe*3)];

    for (int i = 0; i<nbv; i++) {
        head[i] = -1;
    }

    /* Now construct the chain */
    for (int k = 0; k<nbe; k++) {
        for (int j = 0; j<3; j++) {
            int i;
            int p = (3*k)+j;
            i = index[p];
            next[p] = head[i-1];
            head[i-1] = p+1;
            //printf("i=%d, next=%d, head=%d \n", i, next[p], head[i]);
        }
    }

    /* Now we can construct the connectivity maxtrix*/
    int MAXCONNECT = 8;
    int *connectivity = new int[(nbv*MAXCONNECT)];
    int *columns = new int[(nbv*MAXCONNECT)];
    for (int i = 0; i<nbv; i++) {

        /* Go over all of the elements connected to node I*/
        int count = 0;
        int p = head[i];
        while (p != -1) {

            int k = trunc(p/3); //"row in index"
            int j = fmod(p, 3) - 1; //"column" in index

            if (j == -1) {
                j = 2;
                k = k-1;
            }

            /* Sanity check*/
            if (index[((3*k)+j)] != i+1) std::cout << "ERROR - CONNECTIVITY MATRIX" << std::endl;;

            //printf("%d %d %d\n", k, j, i);
            /* Enter element in connectivity matrix */
            connectivity[i*MAXCONNECT + count] = k;
            columns[i*MAXCONNECT + count] = j;
            count++;
            p = next[p-1];
        }
    }
    /*Main loop, allocate a few vectors needed for the computation*/
    double *eta_nbv = new double[nbv];
    double *KVx = new double[nbv];
    double *KVy = new double[nbv];
    double *Eta_nbv = new double[nbe];
    double *kvx = new double[nbe * 3];
    double *kvy = new double[nbe * 3];
    double *dVxdt_abs = new double[nbv];
    double *dVydt_abs = new double[nbv];
    int iter;
    double iterror;
    //Vincent's Section
    //Main loop, allocate a few vectors needed for the computation
    for (iter = 1; iter <= niter; iter++) {

        /*Strain rates - GPU KERNEL 2*/
        derive_xy_elem(dvxdx, dvxdy, vx, index, alpha, beta, nbe);
        derive_xy_elem(dvydx, dvydy, vy, index, alpha, beta, nbe);

        /*'KV' term in equation 22*/
        for (int i = 0; i < nbv; i++) {
            KVx[i] = 0.;
            KVy[i] = 0.;
        }
        for (int n = 0; n < nbe; n++) {

            /*Skip if no ice*/
            if (!isice[n]) continue;

            /*Viscous Deformation*/
            double eta_e = etan[n];
            double eps_xx = dvxdx[n];
            double eps_yy = dvydy[n];
            double eps_xy = .5 * (dvxdy[n] + dvydx[n]);
            for (int i = 0; i < 3; i++) {
                kvx[n * 3 + i] = 2 * Helem[n] * eta_e * (2 * eps_xx + eps_yy) * alpha[n * 3 + i] * areas[n] +
                                 2 * Helem[n] * eta_e * eps_xy * beta[n * 3 + i] * areas[n];
                kvy[n * 3 + i] = 2 * Helem[n] * eta_e * eps_xy * alpha[n * 3 + i] * areas[n] +
                                 2 * Helem[n] * eta_e * (2 * eps_yy + eps_xx) * beta[n * 3 + i] * areas[n];
            }

            /*Add basal friction*/
            if (groundedratio[n] > 0.) {
                int n3 = n * 3;
                double gr_a = groundedratio[n] * areas[n];

                for (int k = 0; k < 3; k++) {
                    for (int i = 0; i < 3; i++) {
                        int i_index = index[n3 + i] - 1;
                        double gr_a_alpha2 = gr_a * alpha2[i_index];

                        for (int j = 0; j < 3; j++) {
                            int j_index = index[n3 + j] - 1;
                            double gr_a_alpha2_vx = gr_a_alpha2 * vx[j_index];
                            double gr_a_alpha2_vy = gr_a_alpha2 * vy[j_index];

                            //     printf("%d, %f, %f, %f \n", n, gr_a, gr_a_alpha2, gr_a_alpha2_vx);
                            if (i == j && j == k) {
                                kvx[n3 + k] = kvx[n3 + k] + gr_a_alpha2_vx / 10.;
                                kvy[n3 + k] = kvy[n3 + k] + gr_a_alpha2_vy / 10.;
                            } else if ((i != j) && (j != k) && (k != i)) {
                                kvx[n3 + k] = kvx[n3 + k] + gr_a_alpha2_vx / 60.;
                                kvy[n3 + k] = kvy[n3 + k] + gr_a_alpha2_vy / 60.;
                            } else {
                                kvx[n3 + k] = kvx[n3 + k] + gr_a_alpha2_vx / 30.;
                                kvy[n3 + k] = kvy[n3 + k] + gr_a_alpha2_vy / 30.;
                            }

                        }

                    }
                }
            }
        }

        for (int i = 0; i < nbv; i++) {
            for (int j = 0; j < 8; j++) {
                if (connectivity[(i * 8 + j)] != 0) {
                    KVx[i] = KVx[i] + kvx[((connectivity[(i * 8 + j)]) - 1) * 3 + ((columns[(i * 8 + j)]))];
                    KVy[i] = KVy[i] + kvy[((connectivity[(i * 8 + j)]) - 1) * 3 + ((columns[(i * 8 + j)]))];
                }
            }
        }
        
        /*Get current viscosity on nodes (Needed for time stepping)*/
        // elem2node(eta_nbv,etan,index,areas,weights,nbe,nbv);

        for (int i = 1; i <= nbe; i++) {
            Eta_nbv[i] = etan[i] * areas[i];
        }


        for (int i = 0; i < nbv; i++) {
            for (int j = 0; j < 8; j++) {
                if (connectivity[(i * 8 + j)] != 0) {
                    eta_nbv[i] = eta_nbv[i] + Eta_nbv[connectivity[(i * 8 + j)] - 1];
                }
            }
        }

        for (int i = 0; i < nbv; i++) { eta_nbv[i] = eta_nbv[i] / weights[i]; }

        /*Velocity rate update in the x and y, refer to equation 19 in Rass paper*/
        double normX = 0.;
        double normY = 0.;

        for (int i = 0; i < nbv; i++) {

            /*1. Get time derivative based on residual (dV/dt)*/
            double ResVx = 1. / (rho * max(80.0, H[i])  * ML[i]) * (-KVx[i] + Fvx[i]); //rate of velocity in the x, equation 23
            double ResVy = 1. / (rho * max(80.0, H[i])  * ML[i]) * (-KVy[i] + Fvy[i]); //rate of velocity in the y, equation 24

            dVxdt[i] = dVxdt[i] * damp + ResVx;
            dVydt[i] = dVydt[i] * damp + ResVy;

            if (isnan(dVxdt[i])) {
                std::cerr << "Found NaN in dVxdt[i]";
                return 1;
            }
            if (isnan(dVydt[i])) {
                std::cerr << "Found NaN in dVydt[i]";
                return 1;
            }

            /*2. Explicit CFL time step for viscous flow, x and y directions*/
            double dtVx = rho * pow(resolx[i], 2) / (4 * eta_nbv[i] * (1. + eta_b) * 4.1);
            double dtVy = rho * pow(resoly[i], 2) / (4 * eta_nbv[i] * (1. + eta_b) * 4.1);


            /*3. velocity update, vx(new) = vx(old) + change in vx, Similarly for vy*/
            vx[i] = vx[i] + relaxation * dVxdt[i] * dtVx;
            vy[i] = vy[i] + relaxation * dVydt[i] * dtVy;

            /*Apply Dirichlet boundary condition,*Residual should also be 0 (for convergence)*/
            if (!isnan(spcvx[i])) {
                vx[i] = spcvx[i];
                dVxdt[i] = 0.;
            }
            if (!isnan(spcvy[i])) {
                vy[i] = spcvy[i];
                dVydt[i] = 0.;
            }

        }


        /*LAST: Update viscosity*/
        for (int i = 0; i < nbe; i++) {
            if (!isice[i]) continue;
            double eps_xx = dvxdx[i];
            double eps_yy = dvydy[i];
            double eps_xy = .5 * (dvxdy[i] + dvydx[i]);
            double EII2 = pow(eps_xx, 2) + pow(eps_yy, 2) + pow(eps_xy, 2) + eps_xx * eps_yy;
            double eta_it = 1.e+14 / 2.;
            if (EII2 > 0.) eta_it = rheology_B[i] / (2 * pow(EII2, (n_glen - 1.) / (2 * n_glen)));

            etan[i] = min(exp(rele * log(eta_it) + (1 - rele) * log(etan[i])), eta_0 * 1e5);
            if (isnan(etan[i])) {
                std::cerr << "Found NaN in etan[i]";
                return 1;
            }
        }

        if ((iter % nout_iter) == 1) {
            double dVxdt_max=0.0, dVydt_max=0.0;
            for (int i = 0; i < nbv; i++){
                dVxdt_max = max(abs(dVxdt[i]), dVxdt_max);
                dVydt_max = max(abs(dVydt[i]), dVydt_max);
            }
            iterror = max(dVxdt_max, dVydt_max);

            std::cout<<"iter="<<iter<<", err="<<iterror<<std::endl;
            if ((iterror < epsi) && (iter > 100)) break;
        }
    }

    cout<< "Iterations: " << iter << "Error: " << iterror << " --> converged" << endl;


    /*Cleanup intermediary vectors*/
    //delete[] eta_nbv;
    //delete[] Eta_nbv;
    delete[] dvxdx;
    delete[] dvxdy;
    delete[] dvydx;
    delete[] dvydy;
    //delete[] KVx;
    //delete[] KVy;
    //delete[] kvx;
    //delete[] kvy;


    /*Write output*/
    fid = fopen(outputfile, "wb");
    if (fid == NULL) std::cerr << "could not open file " << outputfile << " for binary reading or writing";
    WriteData(fid, "PTsolution", "SolutionType");
    WriteData(fid, vx, nbv, 1, "Vx");
    WriteData(fid, vy, nbv, 1, "Vy");
    if (fclose(fid) != 0) std::cerr << "could not close file " << outputfile;



    /*Cleanup and return*/
    delete[] index;
    delete[] x;
    delete[] y;
    delete[] H;
    delete[] surface;
    delete[] base;
    delete[] spcvx;
    delete[] spcvy;
    delete[] ice_levelset;
    delete[] ocean_levelset;
    //delete[] groundedratio;   --munmap_chunk(): invalid pointer error when uncommented (Kevin)
    delete[] isice;
    delete[] rheology_B;
    delete[] rheology_B_temp;
    delete[] vx;
    delete[] vy;
    delete[] friction;
    delete[] alpha2;
    delete[] etan;
    delete[] dVxdt;
    delete[] dVydt;
    delete[] alpha;
    delete[] beta;
    delete[] areas;
    delete[] weights;
    delete[] resolx;
    delete[] resoly;
    delete[] dsdx;
    delete[] dsdy;
    delete[] Helem;
    //delete[] ML;   --free(): invalid next size (noraml) error when uncommented (Kevin)
    //delete[] Fvx;   -- double free or corruption (out) error when uncommented (Kevin)
    //delete[] Fvy;   -- double free or corruption (out) error when uncommented (Kevin)


    return 0;

}/*}}}*/