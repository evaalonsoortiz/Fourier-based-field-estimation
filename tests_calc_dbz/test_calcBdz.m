%% This script uses the phantoms to test calc_dbz.
% The results are compared to the analytical solution, see Chapter 25 of
% Magnetic Resonance Imaging : Physical Principles and Seauence Design,
% Second Edition, Brown R., Cheng Y.,Haacke E., Thompson M., Venkatsen R.

%clearvars;

phantom = "cylinder"
switch(phantom)
%%  An anisotropic rectangular susceptibility in a "little" volume
    case "rect" 
        dim = [16, 16, 16];
        res = [1, 1, 1];
        susin = 1;
        susout = 0;
        sus = zeros(dim) + susout;
        sus(7:10, 8:9, 6:11) = susin;
        if param == 0
        else
            % Translation toward down-right to make the matrix symmetric 
            sus = padarray(sus, [1 1 1], susout, 'pre');
            sus = sus(1:dim(1), 1:dim(2), 1:dim(3));
        end
        
%% A sphere
    case "sphere"
        dim = [256, 256, 256];
        res = [1, 1, 1]; % volume unit
        susin = -0.72e-6; 
        susout = -0.36e-6; 
        radius = 12 % volume unit
        sus_dist = Spherical(dim , res, radius, [susin susout]);
        sus = sus_dist.volume;
%% A cylinder
    case "cylinder"
        dim_without_buffer = [128, 128, 128];
        dim = 1 * dim_without_buffer; % Multiply the dim_without_buffer by a power of 2
        res = [1, 1, 1]; % volume unit
        susin = -0.72e-6; 
        susout = -0.36e-6; 
        radius = 12 % volume unit
        theta = 0 % rad, tilt of the cylinder between B0 and y
        phi = pi/2 % angle between x and measure axis in the xy plane (pi/2 for measure along y, 0 for measure along z)
        sus_dist = Cylindrical(dim, res, radius, theta, [susin susout]);
        sus = sus_dist.volume;

%% A sphere with a bigger volume (add a buffer)
    case "sphere_buffer"
        dim_without_buffer = [128, 128, 128];
        dim = 1*dim_without_buffer; % Multiply by a power of 2
        res = [1, 1, 1]; % volume unit
        susin = -0.72e-6; 
        susout = -0.36e-6;
        radius = 48 % volume unit
        sus_dist = Spherical(dim_without_buffer , res, radius, [susin susout]);
        sus = sus_dist.volume;
        sus = padarray(sus, (dim - dim_without_buffer) / 2, susout, 'post');
        sus = padarray(sus, (dim - dim_without_buffer) / 2, susout, 'pre');
end

% Other parameters
b0 = 1; %[T]
sectionz = round(dim(3) / 2) + 1;
sectiony = round(dim(2) / 2) + 1;
sectionx = round(dim(1) / 2) + 1;

%% Analytical solution

[x,y,z] = ndgrid(linspace(-dim(1)/2, dim(1) / 2, dim(1)), linspace(-dim(2)/2, dim(2) / 2 , dim(2)), linspace(-dim(3)/2, dim(3) / 2, dim(3)));
r = sqrt(x.^2 + y.^2 + z.^2);
tic
if (strcmp(phantom, 'sphere') || strcmp(phantom, 'sphere_buffer'))
    % with Lorentz correction  (p. 753)
    dbz_out = (susin - susout) / 3 .* (radius ./ r).^3 .* (3 .* z.^2 ./ r.^2 - 1) + susout * b0 / 3;
    dbz_out(isnan(dbz_out)) = susout * b0 / 3;
    %dbz_in  = zeros(dim) - susout * b0 / 3;
    dbz_in  = zeros(dim) + (susout + susout*susin) * b0 / (3 + 2 * susout + susin); % Shift from article FBM

    % Create a mask to use the expressions of inside and outside
    % effectively in the sphere an outside
    mask = Spherical(dim, res, radius, [1 0]);
    dbz_in = dbz_in .* mask.volume;
    dbz_out = dbz_out .* (1 - mask.volume);

    dbz_analytical = dbz_in + dbz_out;
    dbz_analytical_ppm = dbz_analytical * 1e6;
        
elseif (strcmp(phantom, 'cylinder'))
    % with Lorentz correction  (p. 753)
    dbz_out = (susin - susout) / 2 .* (radius ./ r).^2 * sin(theta)^2 * cos(2*phi) *b0 + susout * b0 / 3;
    dbz_out(isnan(dbz_out)) = susout * b0 / 3;
    dbz_in  = zeros(dim) + (susin - susout) /6 * b0 * ( 3*cos(theta) - 1)  + susout * b0 / 3;

    % Equivalent to a cylindrical mask because all the measures are
    % done through the center (the axes in the analytical solution and
    % simulation are not identically defined so the cylindrical mask
    % does not suit), EXCEPT FOR THE AXIS PARALLEL TO THE CYLINDER AXES
    mask = Cylindrical(dim, res, radius, theta, [susin susout]);
    dbz_in = dbz_in .* mask.volume;
    dbz_out = dbz_out .* (1 - mask.volume);

    dbz_analytical = dbz_in + dbz_out;
    dbz_analytical_ppm = dbz_analytical * 1e6;
    
end
toc
        
%% Variation calculation
% tic
% dBz_obj = FBFest(sus, res, dim, b0);
% toc
%dBz_map_ppm = real(dBz_obj.volume * 1e6); %TODO remove real ? Loss of the y-translation

tic
    %%---------------------------------------------------------------------- %%
    %% Define constants
    %%---------------------------------------------------------------------- %%

    % k-space window
    k_max = 1./(2.*res);
    interval = 2 * k_max ./ dim;

    % define k-space grid
    [kx,ky,kz] = ndgrid(-k_max(1):interval(1):k_max(1) - interval(1),-k_max(2):interval(2):k_max(2) - interval(2),-k_max(3):interval(3):k_max(3) - interval(3));

    %%---------------------------------------------------------------------- %%
    %% Compute Bdz
    %%---------------------------------------------------------------------- %%

    % compute the fourier transform of the susceptibility distribution
    FT_chi = fftshift(fftn(fftshift(sus)));

    % calculate the scaling coefficient 'kz^2/k^2'
    k2 = kx.^2 + ky.^2 + kz.^2;
    k_scaling_coeff = kz.^2./k2;
    k_scaling_coeff(k2 == 0) = 0;

    % compute Bdz (the z-component of the magnetic field due to a
    % sphere, relative to B0) FFT. The region of interest is
    % assumed to be surrounded by a region of infinite extent whose
    % susceptibility is equal to the susceptibility on the origin
    % corner of the matrix.
    bdzFFT = b0 * (1/3 - k_scaling_coeff).*FT_chi;
    bdzFFT(k2 == 0) = b0 * sus(1, 1, 1) * prod(dim) / 3;
    dbz_volume = ifftshift(ifftn(ifftshift(bdzFFT)));
toc

%dBz_map_ppm = real(volume * 1e6); %TODO remove real ? Loss of the y-translation
dBz_map_ppm = dbz_volume * 1e6; %TODO remove real ? Loss of the y-translation

%% Plot intermediate results
% figure; imagesc(squeeze((1/3-k_scaling_coeff(:, :, sectionz)))); colorbar;
% title('Multiplying volume in FT domain, x section')

k_ifft = ifftshift(ifftn(ifftshift((1/3 - k_scaling_coeff))));

figure; imagesc(squeeze(k_ifft(sectionx, :, :))), colorbar;
title('Ifft of the multiplying volume (spatial domain), x section')

%% Plots to compare with analytical
figure;
plot(squeeze(dbz_analytical_ppm(sectionx, :, sectionz)));
hold on
plot(squeeze(dBz_map_ppm(sectionx, :, sectionz)));
hold off
xlabel('grid position')
ylabel('dBz (ppm)')
legend('Analytical', 'simulation');
title(sprintf('Field in the %s phantom in ppm along y axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))

figure;
plot(squeeze(dbz_analytical_ppm(:, sectiony, sectionz)));
hold on
plot(squeeze(dBz_map_ppm(:, sectiony, sectionz)));
hold off
xlabel('grid position')
ylabel('dBz (ppm)')
legend('Analytical', 'simulation');
title(sprintf('Field in the %s phantom in ppm along x axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))

figure;
plot(squeeze(dbz_analytical_ppm(sectionx, sectiony, :)));
hold on
plot(squeeze(dBz_map_ppm(sectionx, sectiony, :)));
hold off
xlabel('grid position')
ylabel('dBz (ppm)')
legend('Analytical', 'simulation');
title(sprintf('Field in the %s phantom in ppm along z axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))

%% Plot result at a section
%   y section
figure;
subplot(2, 1, 1);
imagesc(squeeze(sus(:, sectiony, :))); colorbar; % colormap winter;
title('susceptibility');

subplot(2, 1, 2);
imagesc(squeeze(dBz_map_ppm(:, sectiony, :))); colorbar;
title('Simulation of the field variation');

sgtitle(sprintf('y section, index %u, %s, radius %u', sectiony, 'sphere', radius))

%% Plots 

% figure;
% plot(squeeze(dbz_analytical_ppm(512/2+1, (512-128)/2+1:(512-128)/2+128, 512/2+1 )), 'LineWidth', 1);
% hold on
% plot(1:128, squeeze(dBz_map_ppm128(128/2+1, :, 128/2+1)), 'LineWidth', 1);
% hold on
% plot(1:128, squeeze(dBz_map_ppm256(256/2+1, (256-128)/2+1:(256-128)/2+128, 256/2+1 )), 'LineWidth', 1);
% hold on
% plot(1:128, squeeze(dBz_map_ppm512(512/2+1, (512-128)/2+1:(512-128)/2+128, 512/2+1 )), 'LineWidth', 1);
% hold off
% xlabel('grid position')
% ylabel('dBz (ppm)')
% legend('Analytical', 'simulation 128^3', 'simulation 256^3', 'simulation 512^3');
% title(sprintf('Field in the %s phantom in ppm along y axis with susin=%0.2e and susout=%0.2e', 'sphere', susin, susout))
% % 
% figure;
% plot(linspace(-256, 256, 512), squeeze(dBz_analytical_ppm512(512/2+1, 512/2+1, :)), 'LineWidth', 1);
% hold on
% plot(linspace(-64, 64, 128), squeeze(dBz_map_ppm128(128/2+1, 128/2+1, :)), 'LineWidth', 1);
% hold on
% plot(linspace(-128, 128, 256), squeeze(dBz_map_ppm256(256/2+1, 256/2+1, :)), 'LineWidth', 1);
% hold on
% plot(linspace(-256, 256, 512), squeeze(dBz_map_ppm512(512/2+1, 512/2+1, :)), 'LineWidth', 1);
% hold off
% xlabel('grid position')
% ylabel('dBz (ppm)')
% legend('Analytical', 'simulation 128^3', 'simulation 256^3', 'simulation 512^3');
% title(sprintf('Field in the %s phantom, radius %u in ppm along z axis with susin=%0.2e and susout=%0.2e', 'sphere', 48, susin, susout))
%
% figure;
% subplot(2, 1, 1)
% plot(squeeze(dbz_analytical_ppm(sectionx, sectiony, :)));
% hold on
% plot(squeeze(dBz_map_ppm(sectionx, sectiony, :)));
% hold off
% xlabel('grid position')
% ylabel('dBz (ppm)')
% legend('Analytical', 'simulation radius 12 dim 128');
% title('Simulation sphere radius 12 dim 128')
% subplot(2, 1, 2)
% plot(squeeze(dBz_map_ppm512(257, 257, :)));
% xlabel('grid position')
% ylabel('dBz (ppm)')
% title('Simulation sphere radius 48 dim 512')
% sgtitle(sprintf('Field in the %s phantom in ppm along z axis with susin=%0.2e and susout=%0.2e', 'sphere', susin, susout))
% 
% figure;
% plot(squeeze(dBz_analytical_ppm512(257, 257, :)), '.-');
% hold on
% plot(linspace(1, 512, 128), squeeze(dBz_map_ppm(sectionx, sectiony, :)), '.-');
% hold on
% plot(squeeze(dBz_map_ppm512(257, 257, :)), '.-');
% hold off
% xlabel('grid position')
% ylabel('dBz (ppm)')
% legend('Analytical', 'simulation radius 12 dim 128','simulation radius 48 dim 512');
% title(sprintf('Field in the %s phantom in ppm along z axis with susin=%0.2e and susout=%0.2e', 'sphere', susin, susout))

% addition = zeros(3*256, 1);
% addition(1:512, 1) = squeeze(dBz_map_ppmr48d512so0(512/2+1, 512/2+1, :));
% addition(129:128+512, 1) = addition(129:128+512, 1) + squeeze(dBz_map_ppm512(512/2+1, 512/2+1, :));
% addition(257:256+512, 1) = addition(257:256+512, 1) + squeeze(dBz_map_ppmr48d512so0(512/2+1, 512/2+1, :));
% 
% figure;
% plot(linspace(-256, 256, 512), squeeze(dBz_map_ppm512(512/2+1, 512/2+1, :)), 'b:', 'LineWidth', 1);
% hold on
% plot(linspace(-128-256, -128+256, 512), squeeze(dBz_map_ppmr48d512so0(512/2+1, 512/2+1, :)), 'b:','LineWidth', 1);
% hold on
% plot(linspace(128-256, 128+256, 512), squeeze(dBz_map_ppmr48d512so0(512/2+1, 512/2+1, :)), 'b:','LineWidth', 1);
% hold on
% plot(linspace(-128-256, 128+256, 3*256), addition, 'b','LineWidth', 1);
% hold on
% plot(linspace(-64, 64, 128), squeeze(dBz_map_ppm128(65, 65, :)), 'LineWidth', 1);
% hold off
% xlabel('grid position')
% ylabel('dBz (ppm)')
% legend('simulation 512^3', 'simulation 512^3 translatée de -128', 'simulation 512^3 translatée de 128', 'addition des simulations + 2 *0,12', 'simulation 128^3');
% title(sprintf('Field in the %s phantom, radius %u in ppm along z axis with susin=%0.2e and susout=%0.2e', 'sphere', 48, susin, susout))

