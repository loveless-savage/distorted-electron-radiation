%% this script generates a 3D model of the beam shape
global Ep w0 z0 tau;
Ep=sqrt(1e18/2.146e18); % sqrt of intensity

w0=2*pi*2; % beam waist
z0=w0^2/2; % focal length
dz=z0; % axial focal misalignment

LaserField=@CircularEnvelope;
x=-40:1.6:40;
y=-40:1.6:40;
z=-100:4:100;
t=0;
E=zeros(51,51,51);
for a=1:50
for b=1:50
for c=1:50
    Zx=z0+i*(z(c)-dz); % coherence length w/ z phase
    Zy=z0+i*(z(c)+dz);
    field=z0/sqrt(Zx*Zy)*exp(-x(a)^2/(2*Zx)-y(b)^2/(2*Zy));
    E(a,b,c)=real(field);
    % Ey=real(x(a)*y(b)/(2*Z^2)*field);
    % Ez=real(-i*x(a)/Z*field);
    % E(c,b,a)=sqrt(Ex^2+Ey^2+Ez^2);
end
end
end
cla
axis vis3d
rotate3d on
isosurface(x,y,z,E,0.2);
view(3);
camlight;
lighting gouraud;