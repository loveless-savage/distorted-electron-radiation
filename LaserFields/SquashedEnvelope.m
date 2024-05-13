%% elliptical gaussian envelope. Two width variables instead of one
function [Ex,Ey,Ez,Bx,By,Bz]=SquashedEnvelope(x,y,z,t)
global Ep w0 z0 tau;

% intermediate constants
squish=4; % by what factor should major+minor axes differ?
wx=w0/sqrt(squish);
wy=w0*sqrt(squish);
zx=wx^2/2; % two focal lengths
zy=wy^2/2;
Zx=zx+i*z; % coherence length w/ z phase
Zy=zy+i*z;
Rx=z+zx^2./z; % radial distance from focus
Ry=z+zy^2./z;
theta=0; % rotate oval counterclockwise by this much
xr=x*cos(theta)-y*sin(theta);
yr=x*sin(theta)+y*cos(theta);

% cross-sectional gaussian envelope
psi=sqrt(zx./Zx*zy./Zy)*exp(-xr.^2./(2*Zx)-yr.^2./(2*Zy));
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

