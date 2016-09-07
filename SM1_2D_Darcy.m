% SM1_2D_Darcy.m applies spectral method
% (Galerkin with Numerical Integration)
% to solve 2d Darcy equation,

%
%   -\alpha u + grad p= f   in \Omega;
%   div u = 0               in \Omega;
%   u.n = k                 on \Gamma_1;
%   p = p_0                 on \Gamma_2;
%
%   domain = rectangular [-1,1]*[-1,1]
%   \Gamma_1=[-1,1]*{1}+{1}*{-1,1}
%   \Gamma_2=[-1,1]*{-1}+{-1}*{-1,1}
%
LY=LX;
NELY = NELX;
dxe = LX/NELX;
dye = dxe;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element

% Generate the Spectral Element mesh
% The domain is partitioned into elements,
% each element contains a cartesian GLL subgrid
[iglob,x,y]=MeshBox(LX,LY,NELX,NELY,NGLL);
nglob = length(x);

% The global numbering of the elements follows this convention:
%
%      ... ... ... ... NELX*NELY
% ^    ... ... ... ... ...
% | NELX+1 ... ... ... 2*NELX
% |     1   2  ... ... NELX  
% --->

% The local numbering of GLL nodes follows this convention:
% 
%  	(1,NGLL)...	(NGLL,NGLL)
% ^   	...	...	... 
% |	(1,1)	...	(NGLL,1)