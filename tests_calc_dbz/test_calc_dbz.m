%% This script uses the phantoms to test calc_dbz.
% The results are compared to the analytical solution, see Chapter 25 of
% Magnetic Resonance Imaging : Physical Principles and Seauence Design,
% Second Edition, Brown R., Cheng Y.,Haacke E., Thompson M., Venkatsen R.

%clearvars;

%% Choose variables for simulation
phantom = "sphere"; % choose phantom shape: sphere, cylinder, ...
field = "demodulated"; % choose field: offset or demodulated
unit = "ppm"; % choose unit: ppm or Hz

switch(phantom)

%% A sphere
    case "sphere"
        subsample_factor = 2; % if = 2: downsampling of the simulation result to achieve the same resolution as from the scanner
        res_factor = 1; % resolution of the phantom [mm]
        
        dim_without_buffer = [128, 128, 128]; % y x z % dimensions of the simulation volume
        dim = [256 256 256]; % buffered dimensions of the volume, used for the fourier transformation
        
        res = res_factor*[1, 1, 1];

        susin = 0; % susceptibility of the material
        susout = 1e-6; % susceptibility of the medium

        % if you only have the susceptibility difference then define value here
        % in this case you can only calculate the demodulated field, if you
        % have the separate values for susin and susout you can also calculate the field offset

        sus_diff = 9e-6; % OR susin - susout;

        radius = 15; % [mm]

        sus_dist = Spherical(dim_without_buffer, res, radius, [sus_diff 0]); % phantom is created
        sus = sus_dist.volume; % extract volume matrix from phantom
        type = 'spherical';
        
%% A cylinder
    case "cylinder"
        subsample_factor = 2; % if = 2: downsampling of the simulation result to achieve the same resolution as from the scanner
        res_factor = 1; % resolution of the phantom [mm]

        dim_without_buffer = [128, 128, 128]; % y x z % dimensions of the simulation volume
        dim = 2*[129, 129, 129]; % buffered dimensions of the volume, used for the fourier transformation
        res = res_factor * [1, 1, 1]; 

        susin = 3e-6; % susceptibility of the material
        susout = -3e-6; % susceptibility of the medium

        % if you only have the susceptibility difference then define value here
        % in this case you can only calculate the demodulated field, if you
        % have the values for susin and susout you can also calculate the field offset
        sus_diff =  6e-6; % OR susin - susout; 
    
        radius = 15; % [mm]
        theta =  pi/2; % rad, tilt of the cylinder between B0 and z, rotation around y-axis
        phi_x = 0; % angle between x and measure axis in the xy plane
        phi_y = pi/2; % angle between y and measure axis in the xy plane

        sus_dist = Cylindrical(dim_without_buffer, res, radius, theta, [sus_diff 0]); % phantom is created
        sus = sus_dist.volume; % extract volume matrix from phantom
        type = 'cylindrical';

end


% Other parameters
b0 = 1; %[T]
gamma_2pi = 42.5775e6; %[Hz/T]

sectionx = round(dim_without_buffer(2) / 2) + 1; % for plotting 
sectiony = round(dim_without_buffer(1) / 2) + 1;
sectionz = round(dim_without_buffer(3) / 2) + 1;

sectionx_subsample = round(dim_without_buffer(2)/subsample_factor / 2) + 1;
sectiony_subsample = round(dim_without_buffer(1)/subsample_factor / 2) + 1;
sectionz_subsample = round(dim_without_buffer(3)/subsample_factor / 2) + 1;

padDim = dim - dim_without_buffer;

%% Analytical solution
res_analytical = res;

[x,y,z] = ndgrid(...
    linspace(-dim_without_buffer(2)/2*res_analytical(2), dim_without_buffer(2) / 2*res_analytical(2), dim_without_buffer(2)), ...
    linspace(-dim_without_buffer(1)/2*res_analytical(1), dim_without_buffer(1) / 2 *res_analytical(1), dim_without_buffer(1)), ...
    linspace(-dim_without_buffer(3)/2*res_analytical(3), dim_without_buffer(3) / 2*res_analytical(3), dim_without_buffer(3)));


r = sqrt(x.^2 + y.^2 + z.^2);
tic
if (strcmp(phantom, 'sphere'))
    % with Lorentz correction  (p. 753)
    dBz_out = (sus_diff) / 3 .* (radius ./ r).^3 .* (3 .* z.^2 ./ r.^2 - 1);
    dBz_out(isnan(dBz_out)) = 0; % because demodulated field, for field offset we add susout/3
    dBz_in  = 0; 
    %dbz_in  = zeros(dim_without_buffer) + (susout + susout*susin) * b0 / (3 + 2 * susout + susin); % Shift from article FBM

    % Create a mask to use the expressions of inside and outside
    % effectively in the sphere an outside
    mask = Spherical(dim_without_buffer, res_analytical, radius, [1 0]);
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
    
    mask = Cylindrical(dim_without_buffer, res_analytical, radius, theta, [1 0]);
    dBz_in = dBz_in .* mask.volume;
    dBz_out_x = dBz_out_x .* (1 - mask.volume);
    dBz_out_y = dBz_out_y .* (1 - mask.volume);
    
    dBz_analytical_x = dBz_in + dBz_out_x; % this is the analytical field offset ratio along x-axis
    dBz_analytical_y = dBz_in + dBz_out_y; % this is the analytical field offset ratio along y-axis
    
end

toc
        
%% Variation calculation
tic 
dBz_obj = FBFest(type, sus, res, dim_without_buffer, sus(1, 1, 1), dim) ; % ( sus, image_res, matrix, sus_ext, b0, dim_with_buff, varargin )

toc
dBz_simulation = dBz_obj.volume; % this is the simulated field offset ratio


% Demodulated field
if (strcmp(field, 'demodulated'))
analytical_x = dBz_analytical_x;
analytical_y = dBz_analytical_y;
simulation = dBz_simulation;
end

% Field offset
if (strcmp(field, 'offset'))
analytical_x = dBz_analytical_x + susout/3;
analytical_y = dBz_analytical_y + susout/3;
simulation = dBz_simulation + susout/3;
end

% Field in Hz
if (strcmp(unit, 'Hz'))
analytical_x = analytical_x * b0 * gamma_2pi;
analytical_y = analytical_y * b0 * gamma_2pi;
simulation = simulation * b0 * gamma_2pi;
end

% Field in ppm
if (strcmp(unit, 'ppm'))
analytical_x = analytical_x * 1e6;
analytical_y = analytical_y * 1e6;
simulation = simulation * 1e6;
end

% Subsample
if subsample_factor == 2
    simulation = sub_sample_3D(simulation, subsample_factor * [1,1,1]);
end


%% Plots to compare with analytical: 
close all
% only along z-axis
figure;
plot(linspace(-dim_without_buffer(3)/2 * res_analytical(3)+ res_analytical(3)/2,dim_without_buffer(3)/2 * res_analytical(3)- res_analytical(3)/2,dim_without_buffer(3)), ...
    squeeze(analytical_x(sectionx, sectiony, :)), ...
    'LineWidth',3,'Color',[0.55 0.55 0.55],'LineStyle',':');
hold on
plot(linspace(-dim_without_buffer(3)/2*res(3)+ res(3)/2,dim_without_buffer(3)/2*res(3) - res(3)/2,dim_without_buffer(3)/subsample_factor), ...
    squeeze(simulation(sectionx_subsample, sectiony_subsample, :)), ...
    'LineWidth',1.5,'Color','r');

hold off
xlabel('z position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
legend('Analytical', 'Simulation');
title(sprintf('Magnetic field %s [%s] in %s (r=%0.1f mm) along z-axis with sus diff=%0.2f ppm', field, unit, phantom, radius, sus_diff*1e6))
grid on

% only along z-axis
figure;
scatter(linspace(-dim_without_buffer(3)/2 *subsample_factor + subsample_factor/2, dim_without_buffer(3)/2*subsample_factor - subsample_factor/2, dim_without_buffer(3)/subsample_factor), ...
    squeeze(simulation(sectionx_subsample, sectiony_subsample, :)), ...
    'filled','LineWidth',0.5,'Color','r');

hold off
xlabel('z position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
%legend('Analytical', 'Simulation');
title(sprintf('Magnetic field %s [%s] in %s (r=%0.1f mm) along z-axis with sus diff=%0.2f ppm', field, unit, phantom, radius, sus_diff*1e6))
grid on

% for different sections together 
figure;
subplot(1, 3, 1);
plot(linspace(-dim_without_buffer(2)/2 * res_analytical(2)+ res_analytical(2)/2, dim_without_buffer(2)/2 * res_analytical(2)- res_analytical(2)/2,dim_without_buffer(2)), ...
    squeeze(analytical_x(sectiony, :, sectionz)), ...
    'LineWidth',3,'Color',[0.55 0.55 0.55],'LineStyle',':');
hold on
plot(linspace(-dim_without_buffer(2)/2*res(2)+ res(3)/2,dim_without_buffer(2)/2*res(2)- res(3)/2,dim_without_buffer(2)/subsample_factor), ...
    squeeze(simulation(sectiony_subsample, :, sectionz_subsample)), ...
    'LineWidth',1.5,'Color','r');

ylim([-5 5])
hold off
xlabel('x position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
legend('Analytical', 'Simulation');
title('Field along x-axis')
grid on

subplot(1, 3, 2);
plot(linspace(-dim_without_buffer(1)/2 * res_analytical(1)+ res_analytical(1)/2,dim_without_buffer(1)/2 * res_analytical(1)- res_analytical(1)/2,dim_without_buffer(1)), ...
    squeeze(analytical_y(:, sectionx, sectionz)), ...
    'LineWidth',3,'Color',[0.55 0.55 0.55],'LineStyle',':');
hold on
plot(linspace(-dim_without_buffer(1)/2*res(1)+ res(3)/2,dim_without_buffer(1)/2*res(1)- res(3)/2,dim_without_buffer(1)/subsample_factor), ...
    squeeze(simulation(:, sectionx_subsample, sectionz_subsample)), ...
    'LineWidth',1.5,'Color','r');

hold off
xlabel('y position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
legend('Analytical', 'Simulation');
title('Field along y-axis')
grid on

subplot(1, 3, 3);
plot(linspace(-dim_without_buffer(3)/2 * res_analytical(3)+ res_analytical(3)/2,dim_without_buffer(3)/2 * res_analytical(3)- res_analytical(3)/2,dim_without_buffer(3)), ...
    squeeze(analytical_x(sectionx, sectiony, :)), ...
    'LineWidth',3,'Color',[0.55 0.55 0.55],'LineStyle',':');
hold on
plot(linspace(-dim_without_buffer(3)/2*res(3)+ res(3)/2,dim_without_buffer(3)/2*res(3)- res(3)/2,dim_without_buffer(3)/subsample_factor), ...
    squeeze(simulation(sectionx_subsample, sectiony_subsample, :)), ...
    'LineWidth',1.5,'Color','r');

hold off
xlabel('z position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
legend('Analytical', 'Simulation');
title('Field along z-axis')
grid on

sgtitle(sprintf('Magnetic field %s [%s] in %s (r=%0.1f mm) with sus diff=%0.2f ppm', field, unit, phantom, radius, sus_diff*1e6))

%% Plot susceptibility difference and simulated field in x, y and z section
figure;
subplot(2, 3, 1);
imagesc(linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    squeeze(sus(:, sectionx, :))*1e6); colorbar; % colormap winter;
title('x-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('z position [mm]')
ylabel('y position [mm]')
axis square

subplot(2, 3, 2);
imagesc(linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    squeeze(sus(sectiony, :, :))*1e6); colorbar; % colormap winter;
title('y-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('z position [mm]')
ylabel('x position [mm]')
axis square

subplot(2, 3, 3);
imagesc(linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    squeeze(sus(:, :, sectionz))*1e6); colorbar; % colormap winter;
title('z-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('x position [mm]')
ylabel('y position [mm]')
axis square

subplot(2, 3, 4);
imagesc(linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    linspace(-dim_without_buffer(2)/2*res(2)+ res(2)/2,dim_without_buffer(2)/2*res(2)- res(2)/2,10), ...
    squeeze(simulation(:, sectionx_subsample, :))); colorbar;
title('x-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('z position [mm]')
ylabel('x position [mm]')
axis square

subplot(2, 3, 5);
imagesc(linspace(-dim_without_buffer(1)/2*res(1)+ res(1)/2,dim_without_buffer(1)/2*res(1)- res(1)/2,10), ...
    linspace(-dim_without_buffer(1)/2*res(1)+ res(1)/2,dim_without_buffer(1)/2*res(1)- res(1)/2,10), ...
    squeeze(simulation(sectiony_subsample, :, :))); colorbar;
title('y-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('z position [mm]')
ylabel('x position [mm]')
axis square

subplot(2, 3, 6);
imagesc(linspace(-dim_without_buffer(3)/2*res(3)+ res(3)/2,dim_without_buffer(3)/2*res(3)- res(3)/2,10), ...
    linspace(-dim_without_buffer(3)/2*res(3)+ res(3)/2,dim_without_buffer(3)/2*res(3)- res(3)/2,10), ...
    squeeze(simulation(:, :, sectionz_subsample))); colorbar;
title('z-section');
xlim([-radius*4 radius*4])
ylim([-radius*4 radius*4])
xlabel('x position [mm]')
ylabel('y position [mm]')
axis square

sgtitle(sprintf('Relative susceptibility maps [ppm] (first row) and simulated field maps [%s] (second row) for x, y and z sections',unit))
