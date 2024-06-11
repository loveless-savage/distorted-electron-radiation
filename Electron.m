clear; close all;
format long
global Ep w0 z0 tau LaserField;
Ep=sqrt(1e18/2.146e18); % sqrt of intensity
w0=2*pi*3; % beam waist
z0=w0^2/2; % focal length
tau=2*pi*12; % pulse duration in radians

% for all available distortions, see the folder ./LaserFields/
addpath("./LaserFields/")
% which distortion are we using? For no distortion use @CircularEnvelope
LaserField = @Chirp;
plottitlestr = join(['Distortion = ',func2str(LaserField)]);

%% sample size & rendering parameters: these will be changed frequently!
% # of electrons
mmax=100;
% polar plot resolution
p_res=100;
% frequency band-pass boundaries.  1=fundamental, 2=second harmonic.
fmin=1.5; %1.3
fmax=2.7; % 2.2
% detector latitude for 2D polar plot
latitude=pi/2;
% which plots to show? default is just a 2D polar plot @ the equator
showFrequencySpectrum=false;
show3dFarField=true;

%% array initialization
% polar phi (longitude) coordinates for plotting
p=0:2*pi/(p_res-1):2*pi;
% accumulate photon counts w/ respect to phi
Np=p*0;
Nt=p*0;
% ...and accumulate those counts across all electrons
Np_tot=p*0;
Nt_tot=p*0;
% accumulate frequency spectra for optional analysis
if showFrequencySpectrum
    Ap_bp_tot=Spectrum([0 1]); % when called w/ only one argument, Spectrum()
    At_bp_tot=Spectrum([0 1]); % ...returns an empty array of appropriate size
end
% phi coordinates converted to x+y
pcos=cos(p);
psin=sin(p);
% latitude converted to cylindrical coordinates
s_lat=sin(latitude);
c_lat=cos(latitude);

%% simulate multiple electrons
for n_e=1:mmax
    % print electron count progress
    n_e

    % random coordinates in the focal plane
    xi=w0;
    yi=w0;
    while xi^2+yi^2>w0^2/2
        xi=w0*2*(rand-0.5);
        yi=w0*2*(rand-0.5);
    end
    zi=z0*2*(rand-0.5);
    % position electron slightly off-center
    %xi=w0/10;
    %yi=0;
    %zi=z0*0;

    % time inverval to compute electron laser interaction
    ti=zi-2*pi*50;
    tf=zi+2*pi*50;
    % initial electron momentum + velocity
    pxi=0;
    pyi=0;
    pzi=0;
    %gamma=sqrt(1+pzi^2);
    %uz=pzi/gamma;
    %h=(1-uz)/(1+uz);
    %h

    % find the electron trajectory
    [t,x,y,z,ux,uy,uz,ax,ay,az]=Trajectory(ti,tf,xi,yi,zi,pxi,pyi,pzi);
    % show that trajectory
    % figure
    % plot(t/(2*pi),x/(2*pi),'b',t/(2*pi),y/(2*pi),'g',t/(2*pi),z/(2*pi),'r');
    
    % scan around the electron for radiation
    for j=1:p_res % phi (longitude)
        % detector cartesian coordinates in 3D
        xd=s_lat*pcos(j);
        yd=s_lat*psin(j);
        zd=c_lat;
        % find electric field vector + time delay at each detector location
        [Ax,Ay,Az,td]=ScatteredField(t,x,y,z,ux,uy,uz,ax,ay,az,xd,yd,zd);
        % dot products to find longitude (At) + latitude (Ap) polarizations
        At=c_lat*pcos(j)*Ax+c_lat*psin(j)*Ay-s_lat*Az;
        Ap=-psin(j)*Ax+pcos(j)*Ay;
        % break into photon frequencies
        [At_bp,nu]=Spectrum(td,At);
        [Ap_bp,nu]=Spectrum(td,Ap);
        % if we are plotting the frequency spectrum, accumulate from each electron
        if showFrequencySpectrum
            At_bp_tot = At_bp_tot+At_bp;
            Ap_bp_tot = Ap_bp_tot+Ap_bp;
        end
        % band-pass frequencies within range of desired harmonic
        for n_y=1:length(nu)
            if nu(n_y)<fmin||nu(n_y)>fmax
                At_bp(n_y)=0;
                Ap_bp(n_y)=0;
            end
        end
        % inverse fourier transform
        ft=trapz(nu,At_bp)*9.00e-4;
        fp=trapz(nu,Ap_bp)*9.00e-4;
        % find number of photons in frequency band
        Nt(j)=ft/(1.55*(fmax+fmin)/2);
        Np(j)=fp/(1.55*(fmax+fmin)/2);
    end
    % accumulate photon counts
    Nt_tot=Nt_tot+Nt;
    Np_tot=Np_tot+Np;
end

%% polar plot of photon counts as a function of phi
polarplot(p,Nt_tot,'g',p,Np_tot,'b','LineWidth',2);
title(plottitlestr);
legend('\theta polarization','\phi polarization','Location','bestoutside');

%% frequency plot
if showFrequencySpectrum
    nw = floor(length(nu)/4);
    plot(nu(1:nw),At_bp_tot(1:nw),'g',nu(1:nw),Ap_bp_tot(1:nw),'b');
    grid on;
end

%% Plot number of photons on entire sphere.
if show3dFarField
    % latitude resolution
    sph_res=p_res;
    [X,Y,Z]=sphere(sph_res);
    % initialize bins in a 2D array based on spherical coordinates
    Nt_sph=zeros(sph_res+1,sph_res+1);
    Np_sph=Nt_sph;
    % paint across all bins
    for k=1:sph_res+1 % theta (latitude)
        k
        for j=1:sph_res+1 % phi (longitude)
            xd=X(k,j);
            yd=Y(k,j);
            zd=Z(k,j);   
            [Ax,Ay,Az,td]=ScatteredField(t,x,y,z,ux,uy,uz,ax,ay,az,xd,yd,zd);
            c_lat=zd;
            s_lat=sqrt(1-c_lat^2);
            pcos=xd/sqrt(xd^2+yd^2+.0001);
            psin=yd/sqrt(xd^2+yd^2+.0001);
            At=c_lat*pcos*Ax+c_lat*psin*Ay-s_lat*Az;
            Ap=-psin*Ax+pcos*Ay;
            At_bp=At.^2;
            Ap2=Ap.^2;
            [At_bp,nu]=Spectrum(td,At);
            [Ap_bp,nu]=Spectrum(td,Ap);
            % find photons in a certain frequency band
            for n_e=1:length(nu)
                if nu(n_e)<fmin||nu(n_e)>fmax
                    At_bp(n_e)=0;
                    Ap_bp(n_e)=0;
                end
            end
            ft=trapz(nu,At_bp)*9.00e-4;
            fp=trapz(nu,Ap_bp)*9.00e-4;
            Nt_sph(k,j)=ft/(1.55*(fmax+fmin)/2);
            Np_sph(k,j)=fp/(1.55*(fmax+fmin)/2);
        end
    end
    % render bins
    fig3d=figure;
    fig3d.Position = [288 338 960 420];
    % phi polarization
    subplot(1,2,1);
    colormap([0 0 0;0 0 .25;0 0 .5;0 0 .75;0 0 1;0 .25 .75;0 .5 .5;0 .75 .25;0 1 0;.25 .75 0;.5 .5 0;.75 .25 .33;1 0 0;1 .25 .25;1 .5 .5;1 .75 .75;1 1 1]);
    surf(X,Y,Z,Np_sph,'Edgecolor', 'none');
    title('\phi Polarization');
    xlabel('x - axis');
    ylabel('y - axis');
    zlabel('z - axis');
    shading interp;
    axis equal;
    alpha(.5);
    hold on;
    zz= -1:1/10:1;
    [xx,yy,zz] = cylinder(w0*zz./z0);
    zz=(zz-.5)*2;
    c=zz*0;
    surf(xx,yy,zz,c,'Edgecolor','white'); % laser cone
    % theta polarization
    subplot(1,2,2);
    colormap([0 0 0;0 0 .25;0 0 .5;0 0 .75;0 0 1;0 .25 .75;0 .5 .5;0 .75 .25;0 1 0;.25 .75 0;.5 .5 0;.75 .25 .33;1 0 0;1 .25 .25;1 .5 .5;1 .75 .75;1 1 1]);
    surf(X,Y,Z,Nt_sph,'Edgecolor', 'none');
    title('\theta Polarization');
    xlabel('x - axis');
    ylabel('y - axis');
    zlabel('z - axis');
    shading interp;
    axis equal;
    alpha(.5);
    hold on;
    surf(xx,yy,zz,c,'Edgecolor','white'); % laser cone
    hold off;
    sgtitle(plottitlestr); % overall title
end

