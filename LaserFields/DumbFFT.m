%% use rudimentary FFT to treat each subfrequency differently in space
function [Ex,Ey,Ez,Bx,By,Bz]=DumbFFT(x,y,z,t)
global Ep w0 z0 tau;

% frequency dependence
N=100; % number of sub-frequencies in the eventual sum
dw=8/tau/N; % frequency step
w=(-N/2:1:N/2)*dw; % range of all frequencies of interest

% intermediate constants
Z=z0+i*z; % coherence length w/ z phase
rho2=x^2+y^2; % sideways radius
R=z+z0^2/z; % radial distance from focus

% spatial gaussian envelope
psi=z0/Z*exp(-(1+w)*rho2/(2*Z));
% frequency envelope
wenv=exp(-(tau/4*w).^2);
% original plane wave for each frequency
plane=exp(i*(1+w)*z);
% scalar field strength in frequency space
wfield=Ep*psi.*wenv.*plane;

% now, we integrate across our given frequency range
field=sum(wfield.*exp(-i*w*t)*dw)*exp(-i*t);

% vector elements
Ex=real(field);
Ey=real(x*y/(2*Z^2)*field);
Ez=real(-i*x/Z*field);
Bx=1*Ey;
By=1*Ex;
Bz=1*real(-i*y/Z*field);