%% use rudimentary FFT to treat each subfrequency differently in space
function [Ex,Ey,Ez,Bx,By,Bz]=DumbFFT(x,y,z,t)
global Ep w0 z0 tau;

% frequency dependence: to make a time envelope of width tau, we
% need a freqency envelope of width 2/tau, which we can cover using
% frequencies within +/-8x the envelope width.
N=100; % number of sub-frequencies in the eventual sum
dw=8/tau*(2/N); % frequency step
wrange=dw*(-N/2:N/2); % range of sub-frequencies
fullenv = zeros(size(t));

% intermediate values for the spatial envelope
Z = z0+i*z; % coherence length w/ z phase
rho2=x.^2+y.^2; % cylindrical radius coordinate
R=z+z0*z0./z; % radial distance from focus
psi=z0./Z.*exp(-rho2./(2*Z)); % beam waist envelope

% build temporal envelope
for n=1:N
    deltaw = wrange(n);
    % basic envelope shape
    env = tau/sqrt(2)*exp(-(tau*deltaw)^2/4);
    % envelope offset
    denv = exp(i*deltaw*(z+rho2./(2*R)-t));
    % add this frequency's envelope to the running sum
    fullenv = fullenv + env.*denv;
end

% plane wave itself
plane = exp(i*(z-t));
field = Ep/sqrt(2*pi).*fullenv.*psi.*plane*dw;

% vector components: Singh model
Ex = real(field);
Ey = real(x.*y./(2*Z.^2).*field);
Ez = real(-i*x./Z.*field);
Bx = 1*Ey;
By = 1*Ex;
Bz = 1*real(-i*y./Z.*field);