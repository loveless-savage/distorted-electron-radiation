function [Ax,Ay,Az,td]=ScatteredField(t,x,y,z,ux,uy,uz,ax,ay,az,xd,yd,zd);
Ax=x*0;
Ay=Ax;
Az=Ax;
td=t*0;
for n=1:length(t);
  cx=xd-ux(n);
  cy=yd-uy(n);
  cz=zd-uz(n);
  gx=cy*az(n)-ay(n)*cz;
  gy=-cx*az(n)+ax(n)*cz;
  gz=cx*ay(n)-ax(n)*cy;
  D=(1-xd*ux(n)-yd*uy(n)-zd*uz(n))^3;
  Ax(n)=-(yd*gz-gy*zd)/D;
  Ay(n)=-(-xd*gz+gx*zd)/D;
  Az(n)=-(xd*gy-gx*yd)/D;
  td(n)=t(n)-x(n)*xd-y(n)*yd-z(n)*zd;
end;