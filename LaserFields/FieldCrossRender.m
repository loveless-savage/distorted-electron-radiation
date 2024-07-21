%% this script takes a cross section of the beam & compares models
global Ep w0 z0 tau;
Ep=sqrt(1e18/2.146e18); % sqrt of intensity

w0=2*pi*2; % beam waist
z0=w0^2/2; % focal length
dz=z0; % axial focal misalignment

LaserField=@CircularEnvelope;