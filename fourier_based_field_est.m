function Bdz = fourier_based_field_est(FOV,matrix,sus,R)

% *************************************************************************
% fourier_based_field_est
%
% DESCRIPTION: Function computes the magnetic field offset produced by a
% sphere (with a given magnetic susceptibility offset relative to the
% surrounding medium) subject to a uniform external magnetic field B0 that
% is orientated along the z-axis.
%
% REFERENCES: J.P. MARQUES, R. BOWTELL Concepts in Magnetic Resonance Part 
%             B (Magnetic Resonance Engineering), Vol. 25B(1) 65?78 (2005)
%
% INPUTS
% FOV = [FOVx FOVy FOVz] in mm
% matrix = [Nx Ny Nz]
% sus = [chi_i chi_e] -> chi_i = internal susceptibility within the sphere
%                                in ppm
%                     -> chi_e = external susceptibility (of the
%                                surrounding medium) in ppm
% R = radius of the sphere in [mm]
%
% OUPUTS
% Bdz = magnetic field offset due to the sphere (relative to B0)
%
% AUTHOR: Eva Alonso Ortiz, PhD (eva.alonso.ortiz@gmail.com)
% DATE LAST MODIFIED: April 2020
%
%*************************************************************************

% =========================== Header ==================================== %
this_fname = 'fourier_based_field_est';
this_info = sprintf('%-20s : ',this_fname);
fprintf([this_info, 'Current date and time: %s\n'], datestr(now));
% ======================================================================= %


%%---------------------------------------------------------------------- %%
%% Define constants 
%%---------------------------------------------------------------------- %%

% image resolution
image_res = FOV./matrix;

% k-space window
k_max = 1./(2.*image_res);

% define image grid
[x,y,z] = ndgrid(linspace(-(matrix(1)-1)/2,(matrix(1)-1)/2,matrix(1)),linspace(-(matrix(2)-1)/2,(matrix(2)-1)/2,matrix(2)),linspace(-(matrix(3)-1)/2,(matrix(3)-1)/2,matrix(3)));

%[x,y,z] = ndgrid(linspace(-(Nx-1)/2,(Nx-1)/2,Nx),linspace(-(Ny-1)/2,(Ny-1)/2,Ny),linspace(-(Nz-1)/2,(Nz-1)/2,Nz));

% define k-space grid
[kx,ky,kz] = ndgrid(linspace(-k_max(1),k_max(1),matrix(1)),linspace(-k_max(2),k_max(2),matrix(2)),linspace(-k_max(3),k_max(3),matrix(3)));

% radial position (in [mm])
r = sqrt((x.*image_res(1)).^2 + (y.*image_res(2)).^2 + (z.*image_res(3)).^2);

% allocate memory for magnetic susceptibility "chi"
chi = zeros(size(x));

chi(r <= R ) = sus(1);
chi(r > R ) = sus(2);

% display suceptibility distrubution
figure;
slice(image_res(2)*y(104:154,104:154,104:154),image_res(1)*x(104:154,104:154,104:154),image_res(3)*z(104:154,104:154,104:154),1e6*chi(104:154,104:154,104:154),0,0,0);
axis equal vis3d; colorbar;
title('\chi [ppm]')
ylabel('y-position [mm]')
xlabel('x-position [mm]')
zlabel('z-position [mm]')

%%---------------------------------------------------------------------- %%
%% Compute Bdz
%%---------------------------------------------------------------------- %%

% compute the fourier transform of the susceptibility distribution
FT_chi = fftshift(fftn(fftshift(chi)));

% calculate the scaling coefficient 'kz^2/k^2' 
k_scaling_coeff = kz.^2./(kx.^2 + ky.^2 + kz.^2);

% set the scaling coefficient to zero when k = 0
k_scaling_coeff((matrix(1)/2+1/2),(matrix(2)/2+1/2),(matrix(3)/2+1/2)) = 0;%-B0*chi_e/3;

% compute Bdz (the z-component of the magnetic field due to a sphere, relative to B0) in [ppm] 
Bdz = ifftshift(ifftn(ifftshift((1/3-k_scaling_coeff).*FT_chi)));


%%---------------------------------------------------------------------- %%
%% Plot some results
%%---------------------------------------------------------------------- %%

% plot saggital view of Bdz
figure;
imagesc(squeeze(1e6*real(Bdz(129,104:154,104:154))));
colorbar;
title('Bdz [ppm]')

xticklabels = squeeze(z(129,129,104))*image_res(3):5:squeeze(z(129,129,154))*image_res(3);
xticks = linspace(1, size(squeeze(1e6*real(Bdz(129,104:154,104:154))), 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('z-position [mm]')

yticklabels = squeeze(y(129,104,129))*image_res(2):5:squeeze(y(129,154,129))*image_res(2);
yticks = linspace(1, size(squeeze(1e6*real(Bdz(129,104:154,104:154))), 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
ylabel('y-position [mm]')

% plot Bdz along the z-axis
figure;
plot(squeeze(z(129,129,104:154))*image_res(3),1e6*squeeze(real(Bdz(129,129,104:154))),'-.k','Linewidth',2);
xlabel('z-position [mm]')
ylabel('B_{dz} [ppm]')

% plot Bdz along the y-axis
figure;
plot(squeeze(y(129,104:154,129))*image_res(2),1e6*squeeze(real(Bdz(129,104:154,129))),'-.k','Linewidth',2);
xlabel('y-position [mm]')
ylabel('B_{dz} [ppm]')
