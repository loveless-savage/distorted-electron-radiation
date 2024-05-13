%% circular gaussian pulse. The golden child of laser foci.
function [Ex,Ey,Ez,Bx,By,Bz]=CircularEnvelope(x,y,z,t)
global Ep z0 tau;

% intermediate constants
Z=z0+i*z; % coherence length w/ z phase
rho2=x.^2+y.^2; % sideways radius
R=z+z0*z0./z; % radial distance from focus

% cross-sectional gaussian envelope
psi=z0./Z.*exp(-rho2./(2*Z));
% temporal pulse envelope
env=exp(-(t-z-rho2./(2*R)).^2/tau^2);
% overall plane behavior
plane=exp(i*(z-t));

% scalar field strength
field=Ep.*psi.*env.*plane;

% vector elements
Ex=real(field);
Ey=real(x.*y./(2*Z.^2).*field);
Ez=real(-i*x./Z.*field);
Bx=1*Ey;
By=1*Ex;
Bz=1*real(-i*y./Z.*field);

