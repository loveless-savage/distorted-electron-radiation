%% use rudimentary FFT to treat each subfrequency differently in space
function [Ex,Ey,Ez,Bx,By,Bz]=DumbFFT(x,y,z,t)
global Ep w0 tau;

% frequency dependence: to make a time envelope of width sqrt(2)*tau, we
% need a freqency envelope of width sqrt(2)/tau, which we can cover using
% frequencies within +/-5x the envelope width.
N=1000; % number of sub-frequencies in the eventual sum
dw=5*2/tau/N; % frequency step
wrange=dw*(-N/2:N/2); % range of sub-frequencies

% to integrate frequencies we will accumulate field vectors
Ex = zeros(size(z));
Ey = zeros(size(z));
Ez = zeros(size(z));
Bx = zeros(size(z));
By = zeros(size(z));
Bz = zeros(size(z));

% use the Singh model to build a gaussian focus for each frequency
for n=1:N
    deltaw = wrange(n);

    % intermediate constants
    dz0=(1+deltaw)*w0^2/2; % focal length
    Z=dz0+i*z; % coherence length w/ z phase
    rho2=x.^2+y.^2; % sideways radius
    R=z+dz0^2./z; % radial distance from focus

    % spatial gaussian envelope
    psi=z0./Z.*exp(-(1+w)*rho2./(2*Z));
    % frequency envelope;
    wenv = exp(-(tau*deltaw)^2/2);
    % plane wave itself
    plane = exp(i*(1+deltaw)*(z-t)); % TODO: (z-t)?
    % complete Singh model beam for this frequency
    wfield = Ep*psi*wenv*plane*dw;

    % integrate vector components
    Ex = Ex + real(wfield);
    Ey = Ey + real(x.*y./(2*Z.^2).*wfield);
    Ez = Ez + real(-i*x./Z.*wfield);
    Bx = Ey;
    By = Ex;
    Bz = Bz + real(-i*y./Z.*wfield);
end
