%% This script uses the phantoms to test calc_dbz.
% The results are compared to the analytical solution, see Chapter 25 of
% Magnetic Resonance Imaging : Physical Principles and Seauence Design,
% Second Edition, Brown R., Cheng Y.,Haacke E., Thompson M., Venkatsen R.

clearvars;

phantom = "sphere_buffer"
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
        dim = [256, 256, 256];
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
        dim_sans_buffer = [256, 256, 256];
        dim = 4 * dim_sans_buffer;
        res = [1, 1, 1]; % volume unit
        susin = -0.72e-6; 
        susout = -0.36e-6; 
        radius = 86 % volume unit
        sus_dist = Spherical(dim_sans_buffer , res, radius, [susin susout]);
        sus = sus_dist.volume;
        sus = padarray(sus, dim_sans_buffer / 2, susout, 'post');
        sus = padarray(sus, dim_sans_buffer / 2, susout, 'pre');
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

    mask = +imutils.masks.spherical_mask(radius, dim, [res, res, res]);
    dbz_in = dbz_in .* mask;
    dbz_out = dbz_out .* (1 - mask);

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
    % does not suit), EXCEPT FOR THE AXIS PARRALLEL TO THE CYLINDER AXES
    mask = +imutils.masks.spherical_mask(radius, dim, [res, res, res]); 
    dbz_in = dbz_in .* mask;
    dbz_out = dbz_out .* (1 - mask);

    dbz_analytical = dbz_in + dbz_out;
    dbz_analytical_ppm = dbz_analytical * 1e6;
    
end
toc
        
%% Variation calculation
tic
dBz_obj = FBFest(sus_dist.volume, sus_dist.image_res, sus_dist.matrix, b0);
toc
dBz_map_ppm = dBz_obj.volume * 1e6;

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
title('Result of the convolution in k-space (before IFFT)');

sgtitle(sprintf('y section, index %u, %s', sectiony, phantom))
