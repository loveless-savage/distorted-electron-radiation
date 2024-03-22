%% spectral decomposition of scattered photons w/ relevant parameters
function [Inu,nu]=Spectrum(t,A)
% resolution
nmax=2^12;
tmax=t(length(t));
tmin=t(1);
% interpolation step sizes
dt=(tmax-tmin)/(nmax-1);
dnu=2*pi/(dt*(nmax-1));

% declare memory for frequency spectrum
nu=0:dnu:(nmax/2)*dnu;
Inu=nu*0;
% in case empty array was requested for accumulation
if nargin==1
  return
end

teven=tmin:dt:tmax;
Aeven=interp1(t,A,teven,'spline');%sample function with even spacing and with points separated by power of 2.
Anu=fft(Aeven)*dt;
for n=1:nmax/2+1
  Inu(n)=abs(Anu(n)).^2/pi;
end