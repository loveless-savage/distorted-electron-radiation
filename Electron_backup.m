clear; close all;
format long
global Ep z0 tau Et;
Ep=sqrt(1e18/2.146e18);
w0=2*pi*4;
z0=w0^2/2;
tau=2*pi*12;
nu1=1.7; %lower frequency detected.  1=fundamental, 2=second harmonic.
nu2=1.86; %upper frequency detected.
mmax=50;
jmax=50;
p=0:2*pi/(jmax-1):2*pi;
Np=p*0;
Nt=p*0;
cp=cos(p);
sp=sin(p);
st=1;
ct=0;

for m=1:mmax
    m
    xi=w0;
    yi=w0;
    while xi^2+yi^2>w0^2
        xi=w0*2*(rand-0.5);
        yi=w0*2*(rand-0.5);
    end
    zi=z0*2*(rand-0.5);
    % xi=w0/10;
    % yi=0;
    % zi=z0*0;
    %time inverval  to compute electron laser interaction
    ti=zi-2*pi*50;
    tf=zi+2*pi*50;
    phi=90*pi/180; %initial electron direction
    theta=180*pi/180; %initial electron direction
    %initial electron momentum
    pxi=0;
    pyi=0;
    pzi=0;
    gamma=sqrt(1+pzi^2);
    uz= pzi/gamma;
    %h=(1-uz)/(1+uz);
    %h
    %find the electron trajectory
    [t,x,y,z,ux,uy,uz,ax,ay,az]=Trajectory(ti,tf,xi,yi,zi,pxi,pyi,pzi);
    % figure
    % plot(t/(2*pi),x/(2*pi),'b',t/(2*pi),y/(2*pi),'g',t/(2*pi),z/(2*pi),'r');
    %
    for j=1:jmax;
        xd=st*cp(j);
        yd=st*sp(j);
        zd=ct;
        [Ax,Ay,Az,td]=ScatteredField(t,x,y,z,ux,uy,uz,ax,ay,az,xd,yd,zd);
        %
        At=ct*cp(j)*Ax+ct*sp(j)*Ay-st*Az;
        Ap=-sp(j)*Ax+cp(j)*Ay;
        [at2,nu]=Spectrum(td,At);
        [ap2,nu]=Spectrum(td,Ap);
        %zero out frequencies outside of band.
        for n=1:length(nu)
            if nu(n)<nu1|nu(n)>nu2
                at2(n)=0;
                ap2(n)=0;
            end
        end
        ft=trapz(nu,at2)*9.00e-4;
        fp=trapz(nu,ap2)*9.00e-4;
        %Find number of photons in frequency band
        Nt1(j)=ft/(1.55*(nu2+nu1)/2);
        Np1(j)=fp/(1.55*(nu2+nu1)/2);
        %
    end
    %
    Nt=Nt+Nt1;
    Np=Np+Np1;
end

figure
polarplot(p,Nt,'g',p,Np,'b','LineWidth',2)


% Plot number of photons on entire sphere.

num =50;
[X,Y,Z]=sphere(num);
Nt=zeros(num+1,num+1);
Np=Nt;
for j=1:num+1;
  j
  for k=1:num+1;
    xd=X(j,k);
    yd=Y(j,k);
    zd=Z(j,k);   
    [Ax,Ay,Az,td]=ScatteredField(t,x,y,z,ux,uy,uz,ax,ay,az,xd,yd,zd);
    %
    %
    ct=zd;
    st=sqrt(1-ct^2);
    cp=xd/sqrt(xd^2+yd^2+.0001);
    sp=yd/sqrt(xd^2+yd^2+.0001);
    At=ct*cp*Ax+ct*sp*Ay-st*Az;
    Ap=-sp*Ax+cp*Ay;
    At2=At.^2;
    Ap2=Ap.^2;
    [at2,nu]=Spectrum(td,At);
    [ap2,nu]=Spectrum(td,Ap);
    %Find photons in a certain frequency band
    for n=1:length(nu)
      if nu(n)<nu1|nu(n)>nu2
        at2(n)=0;
        ap2(n)=0;
      end
    end
    ft=trapz(nu,at2)*9.00e-4;
    fp=trapz(nu,ap2)*9.00e-4;
    Nt(j,k)=ft/(1.55*(nu2+nu1)/2);
    Np(j,k)=fp/(1.55*(nu2+nu1)/2);
  end;
end;
figure;
colormap([0 0 0;0 0 .25;0 0 .5;0 0 .75;0 0 1;0 .25 .75;0 .5 .5;0 .75 .25;0 1 0;.25 .75 0;.5 .5 0;.75 .25 .33;1 0 0;1 .25 .25;1 .5 .5;1 .75 .75;1 1 1]);
surf(X,Y,Z,Np,'Edgecolor', 'none')
xlabel('x - axis')
ylabel('y - axis')
zlabel('z - axis')
shading interp
axis equal
alpha(.5)
hold on
zz= -1:1/10:1;
[xx,yy,zz] = cylinder(w0*zz./z0);
zz=(zz-.5)*2;
c=zz*0;
surf(xx,yy,zz,c,'Edgecolor','white')
hold off
figure;
colormap([0 0 0;0 0 .25;0 0 .5;0 0 .75;0 0 1;0 .25 .75;0 .5 .5;0 .75 .25;0 1 0;.25 .75 0;.5 .5 0;.75 .25 .33;1 0 0;1 .25 .25;1 .5 .5;1 .75 .75;1 1 1]);
surf(X,Y,Z,Nt,'Edgecolor', 'none')
xlabel('x - axis')
ylabel('y - axis')
zlabel('z - axis')
shading interp
axis equal
alpha(.5)
hold on
surf(xx,yy,zz,c,'Edgecolor','white')
hold off


