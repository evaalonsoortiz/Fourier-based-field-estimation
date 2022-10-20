%% This script uses the phantoms to test calc_dbz.
% The results are compared to the analytical solution, see Chapter 25 of
% Magnetic Resonance Imaging : Physical Principles and Seauence Design,
% Second Edition, Brown R., Cheng Y.,Haacke E., Thompson M., Venkatsen R.

%clearvars;

phantom = "cylinder"; % choose if you want a sphere or cylinder phantom
switch(phantom)

%% A sphere
    case "sphere"
        dim_factor = 2;
        dim_without_buffer = dim_factor * [128 128 128];
        buffer_factor = 1;
        dim = buffer_factor * dim_without_buffer;
        res = [1, 1, 1]; % volume unit
        susin = -0.72e-6; 
        susout = -0.36e-6; 
        radius = 12; % volume unit
        sus_dist = Spherical(dim_without_buffer , res, radius, [susin susout]);
        sus = sus_dist.volume;
        type = 'spherical';

%% Sphere with subsampling
    case "sphere_subsample"
        dim_without_buffer = [256 256 256];
        buffer_factor = 1;
        dim = buffer_factor * dim_without_buffer;
        res = [1, 1, 1]; % volume unit
        susin = -0.72e-6; 
        susout = -0.36e-6; 
        radius = 48; % volume unit
        sus_dist = Spherical(dim_without_buffer , res, radius, [susin susout]);
        sus = sub_sample_3D(sus_dist.volume, [2, 2, 2]);
        dim = size(sus); dim_without_buffer = size(sus);
        type = 'spherical';
 
%% A cylinder
    case "cylinder"
        dim_factor = 2;
        dim_without_buffer = dim_factor * [128, 128, 128];
        buffer_factor = 1;
        dim = buffer_factor * dim_without_buffer; % Multiply the dim_without_buffer by a power of 2
        res = [1, 1, 1]; % volume unit
        susin = -0.36e-6; % 1; % -0.72e-6; 
        susout = -9.24e-6; %0; %-0.36e-6; 
        radius = 59; % volume unit
        theta = pi/2; % rad, tilt of the cylinder between B0 and y
        phi = 0; %pi/2 % angle between x and measure axis in the xy plane (pi/2 for measure along y, 0 for measure along z)
        sus_dist = Cylindrical(dim_without_buffer, res, radius, theta, [susin susout]);
        sus = sus_dist.volume;
        type = 'cylindrical';
end

% Other parameters
b0 = 3; %[T]
sectionz = round(dim_without_buffer(3) / 2) + 1;
sectiony = round(dim_without_buffer(2) / 2) + 1;
sectionx = round(dim_without_buffer(1) / 2) + 1;
padDim = dim - dim_without_buffer;

%% Analytical solution
[x,y,z] = ndgrid(linspace(-dim_without_buffer(1)/2, dim_without_buffer(1) / 2, dim_without_buffer(1)), linspace(-dim_without_buffer(2)/2, dim_without_buffer(2) / 2 , dim_without_buffer(2)), linspace(-dim_without_buffer(3)/2, dim_without_buffer(3) / 2, dim_without_buffer(3)));
r = sqrt(x.^2 + y.^2 + z.^2);
tic
if (strcmp(phantom, 'sphere') || strcmp(phantom, 'sphere_subsample'))
    % with Lorentz correction  (p. 753)
    dbz_out = (susin - susout) / 3 .* (radius ./ r).^3 .* (3 .* z.^2 ./ r.^2 - 1) * b0 + susout * b0 / 3;
    dbz_out(isnan(dbz_out)) = susout * b0 / 3;
    %dbz_in  = zeros(dim) - susout * b0 / 3;

    %dbz_in  = zeros(dim_without_buffer) + (susout + susout*susin) * b0 / (3 + 2 * susout + susin); 
    % Shift from article FBM

    dbz_in = susout * b0 / 3;

    % Create a mask to use the expressions of inside and outside
    % effectively in the sphere an
    % outside
    mask = Spherical(dim_without_buffer, res, radius, [1 0]);
    dbz_in = dbz_in .* mask.volume;
    dbz_out = dbz_out .* (1 - mask.volume);

    dbz_analytical = dbz_in + dbz_out;
    
    %ppm and troncate
    %dbz_analytical_ppm = dbz_analytical( padDim(1) + 1:dim_without_buffer(1) + padDim(1), padDim(2) + 1:dim_without_buffer(2) + padDim(2), padDim(3) + 1:dim_without_buffer(3) + padDim(3)) * 1e6;
    dbz_analytical_ppm = dbz_analytical(1:dim_without_buffer(1), 1:dim_without_buffer(2), 1:dim_without_buffer(3)) * 1e6;
        
elseif (strcmp(phantom, 'cylinder'))
    % with Lorentz correction  (p. 753)
    dbz_out = (susin - susout) / 2 .* (radius ./ r).^2 * sin(theta)^2 * cos(2*phi) *b0 + susout * b0 / 3;
    dbz_out(isnan(dbz_out)) = susout * b0 / 3;
    dbz_in  = zeros(dim_without_buffer) + (susin - susout) /6 * b0 * ( 3*cos(theta) - 1)  + susout * b0 / 3;

    % Equivalent to a cylindrical mask because all the measures are
    % done through the center (the axes in the analytical solution and
    % simulation are not identically defined so the cylindrical mask
    % does not suit), EXCEPT FOR THE AXIS PARALLEL TO THE CYLINDER AXES
    mask = Cylindrical(dim_without_buffer, res, radius, theta, [1 0]);
    dbz_in = dbz_in .* mask.volume;
    dbz_out = dbz_out .* (1 - mask.volume);

    dbz_analytical = dbz_in + dbz_out;
    dbz_analytical_ppm = dbz_analytical * 1e6;
    
end
toc
        
%% Variation calculation
tic 
dBz_obj = FBFest(type, sus, res, dim_without_buffer, sus(1, 1, 1), b0, dim); % ( sus, image_res, matrix, sus_ext, b0, dim_with_buff, varargin )
toc
dBz_map_ppm = dBz_obj.volume * 1e6;

%% Plot intermediate results
% figure; imagesc(squeeze((1/3-k_scaling_coeff(:, :, sectionz)))); colorbar;
% title('Multiplying volume in FT domain, x section')

% k_ifft = ifftshift(ifftn(ifftshift((1/3 - k_scaling_coeff))));
% 
% figure; imagesc(squeeze(k_ifft(sectionx, :, :))), colorbar;
% title('Ifft of the multiplying volume (spatial domain), x section')

%% Plots to compare with analytical
close all % close previous plots

figure;
plot(squeeze(dbz_analytical_ppm(sectionx, :, sectionz)));
hold on
plot(squeeze(dBz_map_ppm(sectionx, :, sectionz)));
hold off
xlabel('grid position')
ylabel('dBz (ppm)')
legend('Analytical', 'simulation');
title(sprintf('Field in the %s phantom in ppm along y axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))
grid on

figure;
plot(squeeze(dbz_analytical_ppm(sectionx, sectiony, :)));
hold on
plot(squeeze(dBz_map_ppm(sectionx, sectiony, :)));
hold off
xlabel('grid position')
ylabel('dBz (ppm)')
legend('Analytical', 'simulation');
title(sprintf('Field in the %s phantom in ppm along z axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))
grid on

%% Plot result at a section
%   y section
figure;
subplot(1, 2, 1);
imagesc(squeeze(sus(:, :,sectionz))); colorbar; % colormap winter;
title('susceptibility');
axis square;

subplot(1, 2, 2);
imagesc(squeeze(dBz_map_ppm(:, :, sectionz))); colorbar;
title('Simulation of the field variation');
axis square;

sgtitle(sprintf('y section, index %u, %s, radius %u', sectiony, 'sphere', radius))
