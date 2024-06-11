%% Spatial Chirp: subfrequencies of time pulse spread out side to side
function [Ex,Ey,Ez,Bx,By,Bz]=Chirp(x,y,z,t)
global Ep w0 tau;

% distortion parameters
N=100; % number of sub-frequencies in the eventual sum
chirp_r=0; % sideways frequency displacement
chirp_th=0; % chirp alignment compared to polarization

% frequency dependence: to make a time envelope of width tau, we
% need a freqency envelope of width 2/tau, which we can cover using
% frequencies within +/-8x the envelope width.
dw=8*2/tau/N; % frequency step
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
    R=z+dz0*dz0./z; % radial distance from focus

    % implement chirp displacement
    cx = x + chirp_r/2*cos(chirp_th);
    cy = y + chirp_r/2*cos(chirp_th);
    rho2=cx.^2+cy.^2; % sideways radius

    % spatial gaussian envelope
    psi=dz0./Z.*exp(-(1+deltaw)*rho2./(2*Z));
    % frequency envelope;
    wenv = tau*exp(-(tau*deltaw)^2/4);
    % envelope offset
    denv = exp(i*deltaw*(z+rho2./(2*R)-t));
    % plane wave itself
    plane = exp(i*(z-t));
    % complete Singh model beam for this frequency
    wfield = Ep/sqrt(pi)*psi.*(wenv/2.*denv).*plane*dw;

    % integrate vector components
    Ex = Ex + real(wfield);
    Ey = Ey + real(cx.*cy./(2*Z.^2).*wfield);
    Ez = Ez + real(-i*cx./Z.*wfield);
    Bx = Ey;
    By = Ex;
    Bz = Bz + real(-i*cy./Z.*wfield);
end
