%% two axes that focus at different z distances.
function [Ex,Ey,Ez,Bx,By,Bz]=Astigmatism(x,y,z,t)
global Ep z0 tau;

% intermediate constants
dz=z0/2; % axial focal misalignment
Zx=z0+i*z; % coherence length w/ z phase
Zy=z0+i*(z-dz);
Rx=z+z0^2./z; % radial distance from focus
Ry=z-dz+z0^2./(z-dz);
theta=0; % rotate oval counterclockwise by this much
xr=x*cos(theta)-y*sin(theta);
yr=x*sin(theta)+y*cos(theta);

% cross-sectional gaussian envelope
psi=z0./sqrt(Zx.*Zy).*exp(-xr.^2./(2*Zx)-yr.^2./(2*Zy));
% temporal pulse envelope
env=exp(-(t-z-xr.^2./(2*Rx)-yr.^2./(2*Ry)).^2/tau^2);
% overall plane behavior
plane=exp(i*(z-t));

% scalar field strength
field=Ep.*psi.*plane.*env;

% vector elements (non-rotated x/y)
Ex=real(field);
Ey=real(x.*y./(2*Zx.*Zy).*field);
Ez=real(-i*x./Zx.*field);
Bx=1*Ey;
By=1*Ex;
Bz=1*real(-i*y./Zy.*field);