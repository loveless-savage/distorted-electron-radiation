%% calculate movement of an electron through the laser beam
function [t,x,y,z,ux,uy,uz,ax,ay,az]=Trajectory(ti,tf,xi,yi,zi,pxi,pyi,pzi)
global LaserField;
% qi(1:3) holds 3D position
% qi(4:6) holds 3D momentum
qi=zeros(6,1);
qi(1)=xi;
qi(2)=yi;
qi(3)=zi;
qi(4)=pxi;
qi(5)=pyi;
qi(6)=pzi;
options=odeset('RelTol',1e-9);
% solve the ODE using Lorentz transforms
[t,q]=ode45(@Lorentz,[ti,tf],qi,options);
% q is now an array such that
% q(:,1:3) holds 3D positions
% q(:,4:6) holds 3D momenta
x=q(:,1);
y=q(:,2);
z=q(:,3);
px=q(:,4);
py=q(:,5);
pz=q(:,6);

% allocate space for velocity arrays
ux=x*0;
uy=ux;
uz=ux;
ax=x*0;
ay=ax;
az=ax;
nmax=length(t);
% use electron position & momentum to find velocity & acceleration
for n=1:nmax
   X=x(n);
   Y=y(n);
   Z=z(n);
   T=t(n);
   [Ex,Ey,Ez,Bx,By,Bz]=LaserField(X,Y,Z,T);
   gamma=sqrt(1+px(n)^2+py(n)^2+pz(n)^2);
   ux(n)=px(n)/gamma;
   uy(n)=py(n)/gamma;
   uz(n)=pz(n)/gamma;
   uDotE=ux(n)*Ex+uy(n)*Ey+uz(n)*Ez;
   ax(n)=-(Ex+uy(n)*Bz-By*uz(n)-ux(n)*uDotE)/gamma;
   ay(n)=-(Ey-ux(n)*Bz+Bx*uz(n)-uy(n)*uDotE)/gamma;
   az(n)=-(Ez+ux(n)*By-Bx*uy(n)-uz(n)*uDotE)/gamma;
end