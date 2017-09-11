clear;
close all;

addpath('./lib');
% addpath('./S2 Sampling Toolbox');

% set(0,'defaultAxesFontSize',18);
% set(0,'defaultAxesFontName','Times');
% set(0,'defaultTextFontSize',18);
% set(0,'defaultTextFontName','Times');

%Graphic parameters
WindowLocation;

%% Parameters
%Sound speed [m/s]
c=340.29;

%% Parameters: Reproduced area
%Size of reproduced area [m]
lenX=2.1;
lenY=2.1;

%Interval of reproduced area [m]
dx=0.015;
dy=0.015;

%Number of samples
Nx=lenX/dx;
Ny=lenY/dy;

%Pisitions
xr = repmat(((0:Nx-1)-Nx/2).*dx,Ny,1);
yr = repmat(((0:Ny-1)-Ny/2)'.*dy,1,Nx);
zr = zeros(Nx,Ny);

%% Parameters: Loudspeaker array
%Number of loudspeakers
Nsp = 512;

%Array radius [m]
r_sp = 1.0;

fname_pos = sprintf('./pos/pos%d.mat',Nsp);
if exist(fname_pos,'file')==0 
    [V,Tri,~,Ue]=ParticleSampleSphere('N',Nsp);
    TR=TriRep(Tri,V); 

    figure(81), h=trimesh(TR); set(h,'EdgeColor','b'), axis equal

    [pos.azi,pos.ele,pos.r] = cart2sph(TR.X(:,1),TR.X(:,2),TR.X(:,3));

    save(fname_pos,'pos');
else
    load(fname_pos,'pos');
end

phi_sp = pos.azi;
theta_sp = pi/2-pos.ele;

[x_sp, y_sp, z_sp] = sph2cart(pos.azi,pos.ele,r_sp);

%% Parameters: Microphone array
%Number of microphones
Nm = 32;

%Array radius [m]
r_m = 0.15;

fname_pos = sprintf('./pos/pos%d.mat',Nm);
if exist(fname_pos,'file')==0 
    [V,Tri,~,Ue]=ParticleSampleSphere('N',Nm);
    TR=TriRep(Tri,V); 

    figure(82), h=trimesh(TR); set(h,'EdgeColor','b'), axis equal

    [pos.azi,pos.ele,pos.r] = cart2sph(TR.X(:,1),TR.X(:,2),TR.X(:,3));

    save(fname_pos,'pos');
else
    load(fname_pos,'pos');
end

phi_m = pos.azi;
theta_m = pi/2-pos.ele;

[x_m, y_m, z_m] = sph2cart(pos.azi,pos.ele,r_m);

%% Sound sources

%Source position
rs = 1.5;
theta_s = pi/2;
phi_s = -pi/2;

xs = rs*sin(theta_s)*cos(phi_s);
ys = rs*sin(theta_s)*sin(phi_s);
zs = rs*cos(theta_s);

%Prior location
rs_pr = 1.5;
theta_s_pr = pi/2;
phi_s_pr = -pi/2;

%Amplitude
amp = 1.0;

%Frequency [Hz]
freq=1000;

%Wave number
k=2*pi*freq/c;

%% Parameters

%Maximum order of spherical harmonics
Nmax = ceil(sqrt(Nm)-1);

size_n = Nmax+1;
size_m = 2*Nmax+1;

n = 0:Nmax;
m = -Nmax:Nmax;

%Maximum order for generating signals
Nmax_g = max(ceil(sqrt(Nsp)-1),ceil(k*r_sp));

size_n_g = Nmax_g+1;
size_m_g = 2*Nmax_g+1;

n_g = 0:Nmax_g;
m_g = -Nmax_g:Nmax_g;

%Calculate spherical harmonic functions for higher order
Nmax_c = max(Nmax,Nmax_g);

size_n_c = Nmax_c+1;
size_m_c = 2*Nmax_c+1;

fprintf('Calculating spherical harmonic functions...\n');
Ynm_mat = zeros(size_n_c, size_m_c, Nsp);

Ynm_mat_sp = zeros(size_n_c, size_m_c, Nsp);
for isp=1:Nsp
     Ynm_mat_sp(:,:,isp) = sph_harm_mat(Nmax_c,theta_sp(isp),phi_sp(isp));
end

if Nm==Nsp
    Ynm_mat_m = Ynm_mat_sp;
else
    Ynm_mat_m = zeros(size_n_c, size_m_c, Nm);
    for im=1:Nm
        Ynm_mat_m(:,:,im) = sph_harm_mat(Nmax_c,theta_m(im),phi_m(im));
    end
end

%% Received signals with microphone array
fprintf('Generating received signals...\n');

%Plane wave
sig = zeros(Nm,1);

spec = zeros(size_n_g,size_m_g);
Ynm_g_ps = sph_harm_mat(Nmax_g,theta_s,phi_s);

for in=1:size_n_g
    %Rigid baffle
    spec(in,:) = -k./sph_besselh_diff(n_g(in),1,k.*r_m)*sph_besselh(n_g(in),1,k.*rs).*conj(Ynm_g_ps(in,:));
end

for im=1:Nm
    sig(im) = sum(sum(spec.*Ynm_mat_m(1:size_n_g,1:size_m_g,im)));
end

%% Higher Order Ambisonics

N_hoa = min(ceil(sqrt(Nm)-1),ceil(k*r_m));
size_n_hoa = N_hoa+1;
size_m_hoa = 2*N_hoa+1;

%Encoding stage
Ynm_vec_m = reshape(Ynm_mat_m(1:size_n_hoa,1:size_m_hoa,:),size_n_hoa*size_m_hoa,Nm);

enc_filter_mat = zeros(size_n_hoa,size_m_hoa,Nm);
enc_filter_mat_inv = zeros(size_n_hoa,size_m_hoa,Nm);

enc_filter_spec = -k*r_m^2*sph_besselh_diff(n(1:size_n_hoa),1,k.*r_m)./sph_besselh(n(1:size_n_hoa),1,k.*r_sp);
for im=1:Nm
    enc_filter_mat(:,:,im) = diag(enc_filter_spec)*conj(Ynm_mat_m(1:size_n_hoa,1:size_m_hoa,im));
end
enc_filter_vec = reshape(enc_filter_mat,size_n_hoa*size_m_hoa,Nm);
hoa_enc_vec = enc_filter_vec * sig;

%Decoding stage
Ynm_vec_sp = reshape(conj(Ynm_mat_sp(1:size_n_hoa,1:size_m_hoa,:)),size_n_hoa*size_m_hoa,Nsp);

[U_hoa,s_hoa,V_hoa]=svd(Ynm_vec_sp);
beta = max(max(abs(s_hoa)))*1e-2;

% fprintf('Reg. param. (HOA): %f\n',beta);

hoa_filter_vec = Ynm_vec_sp'/(Ynm_vec_sp*Ynm_vec_sp'+beta*eye(size_n_hoa*size_m_hoa))*enc_filter_vec;

sig_drv_hoa = hoa_filter_vec * sig;

%% Model-based synthesis for prior 

Nmax_mb = min(ceil(sqrt(Nm)-1),ceil(k*r_sp));%ceil(sqrt(Nsp)-1);

size_n_mb = Nmax_mb+1;
size_m_mb = 2*Nmax_mb+1;

n_mb = 0:Nmax_mb;
m_mb = -Nmax_mb:Nmax_mb;

spec_drv_pr = zeros(size_n_mb, size_m_mb);

Ynm_pr = sph_harm_mat(Nmax_mb,theta_s_pr,phi_s_pr);

for in=1:size_n_mb
    spec_drv_pr(in,:) = sph_besselh(n_mb(in),1,k.*rs).*conj(Ynm_pr(in,:))./(2.*pi.*r_sp.*sph_besselh(n_mb(in),1,k.*r_sp));
end

sig_drv_pr = zeros(Nsp,1);
for isp=1:Nsp
    sig_drv_pr(isp) = sum(sum(spec_drv_pr.*Ynm_mat_sp(1:size_n_mb,1:size_m_mb,isp)));
end

%% Source-location-informed recording and reproduction

%Number of control points
Ncp = Nm;

%Array radius [cm]
r_cp = r_m;

phi_cp = phi_m;
theta_cp = theta_m;

%Modified transfer function
H_sp2cp_mod = zeros(Ncp,Nsp);
for icp=1:Ncp
    for isp=1:Nsp
        spec = zeros(size_n_mb,size_m_mb);
        for in=1:size_n_mb
            spec(in,:) = -(1./(k*r_cp^2))*(sph_besselh(n_mb(in),1,k.*r_sp)/sph_besselh_diff(n_mb(in),1,k.*r_cp)).*conj(Ynm_mat_sp(in,1:size_m_mb,isp));
        end
        H_sp2cp_mod(icp,isp) = sum(sum(spec.*Ynm_mat_m(1:size_n_mb,1:size_m_mb,icp)));
    end
end

cov_pr = sig_drv_pr*sig_drv_pr';

%Normalize
cov_pr = cov_pr./trace(cov_pr)*Nsp;

%Power distribution
g_pr = diag(cov_pr);

R = H_sp2cp_mod*diag(g_pr)*H_sp2cp_mod';
R_cor = H_sp2cp_mod*cov_pr*H_sp2cp_mod';

[U_lc,s_lc,V_lc]=svd(R);
alpha = max(max(abs(s_lc)))*1e-2;

% fprintf('Reg. param. (Proposed): %f\n',alpha);

[U,D]= eig(R);
D_diag = diag(D);

Fmap_pw = diag(g_pr)*H_sp2cp_mod'*U*((D+alpha.*eye(Nm))\U');
Fmap_cor = cov_pr*H_sp2cp_mod'*U*((D+alpha.*eye(Nm))\U');

sig_drv_pw = Fmap_pw*sig;
sig_drv_cor = Fmap_cor*sig;

%% Calculate original and reproduced sound field

sig_drv_vec_m1 = sig_drv_hoa;
sig_drv_vec_m2 = sig_drv_pw;
sig_drv_vec_m3 = sig_drv_cor;%sig_drv_pr;%

x_sp_mat = repmat(x_sp',Nx*Ny,1);
y_sp_mat = repmat(y_sp',Nx*Ny,1);
z_sp_mat = repmat(z_sp',Nx*Ny,1);

xr_vec = reshape(xr,Nx*Ny,1);
yr_vec = reshape(yr,Nx*Ny,1);
zr_vec = reshape(zr,Nx*Ny,1);

xr_mat = repmat(xr_vec,1,Nsp);
yr_mat = repmat(yr_vec,1,Nsp);
zr_mat = repmat(zr_vec,1,Nsp);

H_sp2r = green3d(xr_mat, yr_mat, zr_mat, x_sp_mat, y_sp_mat, z_sp_mat, k);
dist_m1 = reshape(H_sp2r*sig_drv_hoa,Nx,Ny);
dist_m2 = reshape(H_sp2r*sig_drv_pw,Nx,Ny);
dist_m3 = reshape(H_sp2r*sig_drv_cor,Nx,Ny);

dist_i = green3d(xr, yr, zr, xs, ys, zs ,k);

%Normalize amplitude
ix_norm = Nx/2+1;
iy_norm = Ny/2+1;

g_i = 1.0/abs(dist_i(ix_norm,iy_norm));
g_m1 = 1.0/abs(dist_m1(ix_norm,iy_norm));
g_m2 = 1.0/abs(dist_m2(ix_norm,iy_norm));
g_m3 = 1.0/abs(dist_m3(ix_norm,iy_norm));

dist_i = g_i.*dist_i;
dist_m1 = g_m1.*dist_m1;
dist_m2 = g_m2.*dist_m2;
dist_m3 = g_m3.*dist_m3;

%Error distribution
err_m1 = 10*log10(abs(dist_m1-dist_i).^2./abs(dist_i).^2);
err_m2 = 10*log10(abs(dist_m2-dist_i).^2./abs(dist_i).^2);
err_m3 = 10*log10(abs(dist_m3-dist_i).^2./abs(dist_i).^2);

%% Draw figures

tt = linspace(0,2*pi,100);

%Sound pressure distribution
zrange = [-4.0 4.0];

fig(1)=figure(1);
set(fig(1),'Position',wloc5_l1(1,:));
hold on;
imagesc(xr_vec,yr_vec,real(dist_i));
plot(r_sp*cos(tt), r_sp*sin(tt), 'w', 'LineWidth',4);
plot(r_sp*cos(tt), r_sp*sin(tt), 'k', 'LineWidth',2);
hold off;
axis equal;
axis tight;
caxis(zrange);
colormap(flipud(pink));
colorbar;
xlabel('x [m]'); ylabel('y [m]');

fig(2)=figure(2);
set(fig(2),'Position',wloc5_l1(2,:));
hold on;
imagesc(xr_vec,yr_vec,real(dist_m1));
plot(r_sp*cos(tt), r_sp*sin(tt), 'w', 'LineWidth',4);
plot(r_sp*cos(tt), r_sp*sin(tt), 'k', 'LineWidth',2);
hold off;
axis equal;
axis tight;
caxis(zrange);
colormap(flipud(pink));
colorbar;
xlabel('x [m]'); ylabel('y [m]');

fig(3)=figure(3);
set(fig(3),'Position',wloc5_l1(3,:));
hold on;
imagesc(xr_vec,yr_vec,real(dist_m2));
plot(r_sp*cos(tt), r_sp*sin(tt), 'w', 'LineWidth',4);
plot(r_sp*cos(tt), r_sp*sin(tt), 'k', 'LineWidth',2);
hold off;
axis equal;
axis tight;
caxis(zrange);
colormap(flipud(pink));
colorbar;
xlabel('x [m]'); ylabel('y [m]');

fig(4)=figure(4);
set(fig(4),'Position',wloc5_l1(4,:));
hold on;
imagesc(xr_vec,yr_vec,real(dist_m3));
plot(r_sp*cos(tt), r_sp*sin(tt), 'w', 'LineWidth',4);
plot(r_sp*cos(tt), r_sp*sin(tt), 'k', 'LineWidth',2);
hold off;
axis equal;
axis tight;
caxis(zrange);
colormap(flipud(pink));
colorbar;
xlabel('x [m]'); ylabel('y [m]');

%Error distribution
zrange = [-30, 0];
contour_vec = zrange(1):3:zrange(2);

fig(12)=figure(12);
set(fig(12),'Position',wloc5_l2(2,:));
hold on;
contourf(xr,yr,err_m1,contour_vec);
plot(r_sp*cos(tt), r_sp*sin(tt), 'w', 'LineWidth',4);
plot(r_sp*cos(tt), r_sp*sin(tt), 'k', 'LineWidth',2);
hold off;
axis equal;
axis tight;
caxis(zrange);
colormap(flipud(pink));
colorbar;
xlabel('x [m]'); ylabel('y [m]');

fig(13)=figure(13);
set(fig(13),'Position',wloc5_l2(3,:));
hold on;
contourf(xr,yr,err_m2,contour_vec);
plot(r_sp*cos(tt), r_sp*sin(tt), 'w', 'LineWidth',4);
plot(r_sp*cos(tt), r_sp*sin(tt), 'k', 'LineWidth',2);
hold off;
axis equal;
axis tight;
caxis(zrange);
colormap(flipud(pink));
colorbar;
xlabel('x [m]'); ylabel('y [m]');

fig(14)=figure(14);
set(fig(14),'Position',wloc5_l2(4,:));
hold on;
contourf(xr,yr,err_m3,contour_vec);
plot(r_sp*cos(tt), r_sp*sin(tt), 'w', 'LineWidth',4);
plot(r_sp*cos(tt), r_sp*sin(tt), 'k', 'LineWidth',2);
hold off;
axis equal;
axis tight;
caxis(zrange);
colormap(flipud(pink));
colorbar;
xlabel('x [m]'); ylabel('y [m]');