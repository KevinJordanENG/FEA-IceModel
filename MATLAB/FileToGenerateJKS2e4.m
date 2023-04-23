md = loadmodel ('./JKS2e4');
% ======================== LOAD DATA FROM MD ===============================

nbe              = md.mesh.numberofelements;
nbv              = md.mesh.numberofvertices;
g                = md.constants.g;
rho              = md.materials.rho_ice;
rho_w            = md.materials.rho_water;
yts              = md.constants.yts;
index            = md.mesh.elements;
spcvx            = md.stressbalance.spcvx/yts;
spcvy            = md.stressbalance.spcvy/yts;
x                = md.mesh.x;
y                = md.mesh.y;
H                = md.geometry.thickness;
surface          = md.geometry.surface;
base             = md.geometry.base;
ice_levelset     = md.mask.ice_levelset;
ocean_levelset   = md.mask.ocean_levelset;
rheology_B_temp  = md.materials.rheology_B;
vx               = md.initialization.vx/yts;
vy               = md.initialization.vy/yts;
friction         = md.friction.coefficient;

% =================== SIMILAR TO C++ FILE LINE 384 DOWN ===================


	[alpha beta gamma]=GetNodalFunctionsCoeff(index,x,y);
	areas=GetAreas(index,x,y);
 % }}}

save 2DSSA