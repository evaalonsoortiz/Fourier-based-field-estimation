%% This script uses the phantoms to test calc_dbz.

% run this script in the main folder (with FBFest.m) after adding the folder test_calc_dbz to the path

% The results are compared to the analytical solution, see Chapter 25 of
% Magnetic Resonance Imaging : Physical Principles and Seauence Design,
% Second Edition, Brown R., Cheng Y.,Haacke E., Thompson M., Venkatsen R.

%clearvars;

%% Choose variables for simulation
phantom = "sphere"; % choose phantom shape: sphere, cylinder, ...
field = "demodulated"; % choose field: offset or demodulated, default is demodulated
unit = "ppm"; % choose unit: ppm or Hz, default is ppm

% Other parameters
b0 = 1; % B0 field strength [T]
gamma_2pi = 42.5775; % Larmor frequency of Hydrogen atoms [MHz/T]

switch(phantom)

%% A sphere
    case "sphere"
        radius = 15; % [mm]

        % define value of susceptibility difference as specific value or calculate susin (chi_i) - susout (chi_e)
        % if you have separate values then you can also calculate the field offset, otherwise only the demodulated field is possible

        %susin = 0; % susceptibility of the material [ppm]
        %susout = 1; % susceptibility of the medium [ppm]
        sus_diff = 9; % OR susin - susout; susceptibility difference [ppm]

        dim_without_buffer = [128, 128, 128]; % y x z % dimensions of the simulation volume
        dim = 2 * dim_without_buffer; % buffered dimensions of the volume, used for the fourier transformation
        
        res = [1, 1, 1]; % resolution of the phantom [mm]

        subsample_factor = 1; % if = 2: downsampling of the simulation result to achieve the same resolution as from the scanner

        sus_dist = Spherical(dim_without_buffer, res, radius, [sus_diff 0]); % phantom is created
        sus = sus_dist.volume; % extract volume matrix from phantom
        type = 'spherical';
        
%% A cylinder
    case "cylinder"
        radius = 15; % [mm]
        theta =  pi/2; % angle [rad] of rotation of the cylinder axis (starts parallel to z-axis and $B_0$) around y-axis
        phi_x = 0; % angle [rad] between x and measure axis in the xy plane
        phi_y = pi/2; % angle [rad] between y and measure axis in the xy plane

        % define value of susceptibility difference as specific value or calculate susin (chi_i) - susout (chi_e)
        % if you have separate values then you can also calculate the field offset, otherwise only the demodulated field is possible

        %susin = 0; % susceptibility of the material [ppm]
        %susout = 1; % susceptibility of the medium [ppm]
        sus_diff = 9; % OR susin - susout; susceptibility difference [ppm]

        dim_without_buffer = [128, 128, 128]; % y x z % dimensions of the simulation volume
        dim = 2 * dim_without_buffer; % buffered dimensions of the volume, used for the fourier transformation

        res = [1, 1, 1]; % resolution of the phantom [mm]

        subsample_factor = 1; % if = 2: downsampling of the simulation result to achieve the same resolution as from the scanner

        sus_dist = Cylindrical(dim_without_buffer, res, radius, theta, [sus_diff 0]); % phantom is created
        sus = sus_dist.volume; % extract volume matrix from phantom
        type = 'cylindrical';

end


%% Analytical solution

[x,y,z] = ndgrid(...
    linspace(-dim_without_buffer(2)/2*res(2), dim_without_buffer(2) / 2*res(2), dim_without_buffer(2)), ...
    linspace(-dim_without_buffer(1)/2*res(1), dim_without_buffer(1) / 2 *res(1), dim_without_buffer(1)), ...
    linspace(-dim_without_buffer(3)/2*res(3), dim_without_buffer(3) / 2*res(3), dim_without_buffer(3)));


r = sqrt(x.^2 + y.^2 + z.^2);
tic
if (strcmp(phantom, 'sphere'))
    % with Lorentz correction  (p. 753)
    dBz_out = (sus_diff) / 3 .* (radius ./ r).^3 .* (3 .* z.^2 ./ r.^2 - 1);
    dBz_out(isnan(dBz_out)) = 0; % because demodulated field, for field offset we add susout/3
    dBz_in  = 0; 

    % Create a mask to use the expressions of inside and outside
    % effectively in the sphere an outside
    mask = Spherical(dim_without_buffer, res, radius, [1 0]);
    dBz_in = dBz_in .* mask.volume;
    dBz_out = dBz_out .* (1 - mask.volume);

    dBz_analytical_x = dBz_in + dBz_out; % this is the analytical field offset ratio
    dBz_analytical_y = dBz_in + dBz_out; 

        
elseif (strcmp(phantom, 'cylinder'))
    % with Lorentz correction  (p. 753)
    % along x-axis phi = 0
    dBz_out_x = (sus_diff) / 2 .* (radius ./ r).^2 * sin(theta)^2 * cos(2*phi_x);
    dBz_out_x(isnan(dBz_out_x)) = 0;

    % along y-axis phi = pi/2
    dBz_out_y = (sus_diff) / 2 .* (radius ./ r).^2 * sin(theta)^2 * cos(2*phi_y);
    dBz_out_y(isnan(dBz_out_x)) = 0;

    dBz_in  = zeros(dim_without_buffer) + (sus_diff) /6 * ( 3*cos(theta) - 1);

    % Equivalent to a cylindrical mask because all the measures are
    % done through the center (the axes in the analytical solution and
    % simulation are not identically defined so the cylindrical mask
    % does not suit), EXCEPT FOR THE AXIS PARALLEL TO THE CYLINDER AXES
    
    mask = Cylindrical(dim_without_buffer, res, radius, theta, [1 0]);
    dBz_in = dBz_in .* mask.volume;
    dBz_out_x = dBz_out_x .* (1 - mask.volume);
    dBz_out_y = dBz_out_y .* (1 - mask.volume);
    
    dBz_analytical_x = dBz_in + dBz_out_x; % this is the analytical field offset ratio along x-axis
    dBz_analytical_y = dBz_in + dBz_out_y; % this is the analytical field offset ratio along y-axis
    
end

toc
        
%% Variation calculation
tic 
dBz_obj = FBFest(type, sus, res, dim_without_buffer, dim) ;

toc
dBz_simulation = dBz_obj.volume; % this is the simulated field offset ratio

%% Subsample simulation to lower resolution 
% (only symmetrical for factor 2)
dBz_simulation = sub_sample_3D(dBz_simulation, subsample_factor * [1,1,1]);


%% Change to field offset and/or Hz
% Default is demodulated field in ppm

% Field offset
if (strcmp(field, 'offset'))
dBz_analytical_x = dBz_analytical_x + susout/3;
dBz_analytical_y = dBz_analytical_y + susout/3;
dBz_simulation = dBz_simulation + susout/3;
end

% Field in Hz
if (strcmp(unit, 'Hz'))
dBz_analytical_x = dBz_analytical_x * b0 * gamma_2pi;
dBz_analytical_y = dBz_analytical_y * b0 * gamma_2pi;
dBz_simulation = dBz_simulation * b0 * gamma_2pi;
end



%% Plots to compare with analytical: 

sectionx = round(dim_without_buffer(2) / 2) + 1; % for plotting 
sectiony = round(dim_without_buffer(1) / 2) + 1;
sectionz = round(dim_without_buffer(3) / 2) + 1;

sectionx_subsample = round(dim_without_buffer(2)/subsample_factor / 2) + 1;
sectiony_subsample = round(dim_without_buffer(1)/subsample_factor / 2) + 1;
sectionz_subsample = round(dim_without_buffer(3)/subsample_factor / 2) + 1;

max_val = max([max(dBz_analytical_x, [],'all'), max(dBz_analytical_y, [],'all'), max(dBz_simulation, [],'all')]);
min_val = min([min(dBz_analytical_x, [],'all'), min(dBz_analytical_y, [],'all'), min(dBz_simulation, [],'all')]);
max_lim = max_val + 0.1 * max_val;
min_lim = min_val + 0.1 * min_val;

close all

% plot along x, y and z axis of the simulated and analytical field
figure;
subplot(1, 3, 1);
plot(linspace(-size(dBz_analytical_x,2)/2 * res(2)+ res(2)/2, size(dBz_analytical_x,2)/2 * res(2)- res(2)/2,size(dBz_analytical_x,2)), ...
    squeeze(dBz_analytical_x(sectiony, :, sectionz)), ...
    'LineWidth',2,'Color','r');
hold on
plot(linspace(-size(dBz_simulation,2)/2 * subsample_factor + subsample_factor/2,size(dBz_simulation,2)/2*subsample_factor - subsample_factor/2,size(dBz_simulation,2)), ...
    squeeze(dBz_simulation(sectiony_subsample, :, sectionz_subsample)), ...
    'LineWidth',2,'Color','k','LineStyle',':');

ylim([min_lim max_lim])
hold off
xlabel('x position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
legend('Analytical', 'Simulation');
title('Field along x-axis')
grid on

subplot(1, 3, 2);
plot(linspace(-size(dBz_analytical_y,1)/2 * res(1)+ res(1)/2,size(dBz_analytical_y,1)/2 * res(1)- res(1)/2,size(dBz_analytical_y,1)), ...
    squeeze(dBz_analytical_y(:, sectionx, sectionz)), ...
    'LineWidth',2,'Color','r');
hold on
plot(linspace(-size(dBz_simulation,1)/2 * subsample_factor + subsample_factor/2,size(dBz_simulation,1)/2*subsample_factor - subsample_factor/2,size(dBz_simulation,1)), ...
    squeeze(dBz_simulation(:, sectionx_subsample, sectionz_subsample)), ...
    'LineWidth',2,'Color','k','LineStyle',':');

ylim([min_lim max_lim])
hold off
xlabel('y position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
legend('Analytical', 'Simulation');
title('Field along y-axis')
grid on

subplot(1, 3, 3);
plot(linspace(-size(dBz_analytical_x,3)/2 * res(3)+ res(3)/2,size(dBz_analytical_x,3)/2 * res(3)- res(3)/2,size(dBz_analytical_x,3)), ...
    squeeze(dBz_analytical_x(sectionx, sectiony, :)), ...
    'LineWidth',2,'Color','r');
hold on
plot(linspace(-size(dBz_simulation,3)/2 * subsample_factor + subsample_factor/2,size(dBz_simulation,3)/2*subsample_factor - subsample_factor/2,size(dBz_simulation,3)), ...
    squeeze(dBz_simulation(sectiony_subsample, sectionx_subsample, :)), ...
    'LineWidth',2,'Color','k','LineStyle',':');

ylim([min_lim max_lim])
hold off
xlabel('z position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
legend('Analytical', 'Simulation');
title('Field along z-axis')
grid on

sgtitle(sprintf('Magnetic field %s [%s] in %s (r=%0.1f mm) with sus diff=%0.2f ppm', field, unit, phantom, radius, sus_diff))


%% Plot susceptibility difference and simulated field in x, y and z section
figure;
subplot(2, 3, 1);
imagesc(linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    squeeze(sus(:, sectionx, :))); colorbar; % colormap winter;
title('x-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('z position [mm]')
ylabel('y position [mm]')
axis square

subplot(2, 3, 2);
imagesc(linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    squeeze(sus(sectiony, :, :))); colorbar; % colormap winter;
title('y-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('z position [mm]')
ylabel('x position [mm]')
axis square

subplot(2, 3, 3);
imagesc(linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    squeeze(sus(:, :, sectionz))); colorbar; % colormap winter;
title('z-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('x position [mm]')
ylabel('y position [mm]')
axis square

subplot(2, 3, 4);
imagesc(linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    squeeze(dBz_simulation(:, sectionx_subsample, :)), [min_lim, max_lim]); colorbar;
title('x-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('z position [mm]')
ylabel('y position [mm]')
axis square

subplot(2, 3, 5);
imagesc(linspace(-dim_without_buffer(1)/2*res(1)+ res(1)/2,dim_without_buffer(1)/2*res(1)- res(1)/2,10), ...
    linspace(-dim_without_buffer(1)/2*res(1)+ res(1)/2,dim_without_buffer(1)/2*res(1)- res(1)/2,10), ...
    squeeze(dBz_simulation(sectiony_subsample, :, :)), [min_lim, max_lim]); colorbar;
title('y-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('z position [mm]')
ylabel('x position [mm]')
axis square

subplot(2, 3, 6);
imagesc(linspace(-dim_without_buffer(3)/2*res(3)+ res(3)/2,dim_without_buffer(3)/2*res(3)- res(3)/2,10), ...
    linspace(-dim_without_buffer(3)/2*res(3)+ res(3)/2,dim_without_buffer(3)/2*res(3)- res(3)/2,10), ...
    squeeze(dBz_simulation(:, :, sectionz_subsample)), [min_lim, max_lim]); colorbar;
title('z-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('x position [mm]')
ylabel('y position [mm]')
axis square

sgtitle(sprintf('Relative susceptibility maps [ppm] (first row) and simulated field maps [%s] (second row) for x, y and z sections',unit))
