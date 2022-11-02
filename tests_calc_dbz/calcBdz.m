%% This script uses the phantoms to test calc_dbz.
% The results are compared to the analytical solution, see Chapter 25 of
% Magnetic Resonance Imaging : Physical Principles and Seauence Design,
% Second Edition, Brown R., Cheng Y.,Haacke E., Thompson M., Venkatsen R.

%clearvars;

%% Choose variables for simulation
phantom = "cylinder"; % choose phantom shape: rect, sphere, cylinder, ...
field = "demodulated"; % choose field: offset or demodulated
unit = "Hz"; % choose unit: ppm or Hz



switch(phantom)
%%  An anisotropic rectangular susceptibility in a "little" volume
    case "rect" 
        dim = [16, 16, 16];
        res = [1, 1, 1];
        susin = 1;
        susout = 0;
        sus = zeros(dim) + susout;
        sus(7:10, 8:9, 6:11) = susin;
        type = '';
        
        
%% A sphere
    case "sphere"
        dim = [256, 256, 256];
        dim_without_buffer = dim;
        res = [1, 1, 1]; % volume unit
        susin = 0; %-0.72e-6; 
        susout = 1; %-0.36e-6; 
        radius = 10; % volume unit
        sus_dist = Spherical(dim , res, radius, [susin susout]);
        sus = sus_dist.volume;
        type = 'spherical';
        
%% A cylinder
    case "cylinder"
        res_factor = 1;
        dim_without_buffer = 2*[128, 128, 128];
        dim = [256, 256, 256]; 
        res = res_factor * [1, 1, 1]; % volume unit
        susin = 3.60e-7; % -0.72e-6; 
        susout = -9.24e-6; %-0.36e-6; 
        radius = 5.9; % volume unit
        theta =  pi/2; % rad, tilt of the cylinder between B0 and y
        phi = 0; %pi/2 % angle between x and measure axis in the xy plane (pi/2 for measure along y, 0 for measure along x)
        sus_dist = Cylindrical(dim_without_buffer, res, radius, theta, [susin susout]);
        sus = sus_dist.volume;
        type = 'cylindrical';

end


% Other parameters
b0 = 1; %[T]
gamma_2pi = 42.5775e6; %[Hz/T]
sectionz = round(dim_without_buffer(3) / 2) + 1;
sectiony = round(dim_without_buffer(2) / 2) + 1;
sectionx = round(dim_without_buffer(1) / 2) + 1;
padDim = dim - dim_without_buffer;

%% Analytical solution

%[x,y,z] = ndgrid(linspace(-dim_without_buffer(1)/2, dim_without_buffer(1) / 2, dim_without_buffer(1)), linspace(-dim_without_buffer(2)/2, dim_without_buffer(2) / 2 , dim_without_buffer(2)), linspace(-dim_without_buffer(3)/2, dim_without_buffer(3) / 2, dim_without_buffer(3)));
[x,y,z] = ndgrid(linspace(-dim_without_buffer(1)/2 * res(1), dim_without_buffer(1) / 2 * res(1), dim_without_buffer(1)), ...
    linspace(-dim_without_buffer(2)/2 * res(2), dim_without_buffer(2) / 2  * res(2), dim_without_buffer(2)), ...
    linspace(-dim_without_buffer(3)/2 * res(3), dim_without_buffer(3) / 2 * res(3), dim_without_buffer(3)));
r = sqrt(x.^2 + y.^2 + z.^2);
tic
if (strcmp(phantom, 'sphere'))
    % with Lorentz correction  (p. 753)
    dBz_out = (susin - susout) / 3 .* (radius ./ r).^3 .* (3 .* z.^2 ./ r.^2 - 1) + susout / 3;
    dBz_out(isnan(dBz_out)) = susout / 3;
    dBz_in  = susout / 3;
    %dbz_in  = zeros(dim_without_buffer) + (susout + susout*susin) * b0 / (3 + 2 * susout + susin); % Shift from article FBM

    % Create a mask to use the expressions of inside and outside
    % effectively in the sphere an outside
    mask = Spherical(dim_without_buffer, res, radius, [1 0]);
    dBz_in = dBz_in .* mask.volume;
    dBz_out = dBz_out .* (1 - mask.volume);

    dBz_analytical = dBz_in + dBz_out; % this is the analytical field offset ratio
    
    %ppm and troncate
    %dbz_analytical_ppm = dbz_analytical( padDim(1) + 1:dim_without_buffer(1) + padDim(1), padDim(2) + 1:dim_without_buffer(2) + padDim(2), padDim(3) + 1:dim_without_buffer(3) + padDim(3)) * 1e6;
    
        
elseif (strcmp(phantom, 'cylinder'))
    % with Lorentz correction  (p. 753)
    dBz_out = (susin - susout) / 2 .* (radius ./ r).^2 * sin(theta)^2 * cos(2*phi) + susout / 3;
    dBz_out(isnan(dBz_out)) = susout / 3;
    dBz_in  = zeros(dim_without_buffer) + (susin - susout) /6 * ( 3*cos(theta) - 1)  + susout / 3;

    % Equivalent to a cylindrical mask because all the measures are
    % done through the center (the axes in the analytical solution and
    % simulation are not identically defined so the cylindrical mask
    % does not suit), EXCEPT FOR THE AXIS PARALLEL TO THE CYLINDER AXES
    mask = Cylindrical(dim_without_buffer, res, radius, theta, [1 0]);
    dBz_in = dBz_in .* mask.volume;
    dBz_out = dBz_out .* (1 - mask.volume);
    
    dBz_analytical = dBz_in + dBz_out; % this is the analytical field offset ratio
    
end

toc
        
%% Variation calculation
tic 
dBz_obj = FBFest(type, sus, res, dim_without_buffer, sus(1, 1, 1), dim) ; % ( sus, image_res, matrix, sus_ext, b0, dim_with_buff, varargin )

toc
dBz_simulation = dBz_obj.volume; % this is the simulated field offset ratio


% Demodulated field
if (strcmp(field, 'demodulated'))
analytical = dBz_analytical - susout/3;
simulation = dBz_simulation - susout/3;
end

% Field offset
if (strcmp(field, 'offset'))
analytical = dBz_analytical;
simulation = dBz_simulation;
end

% Field in Hz
if (strcmp(unit, 'Hz'))
analytical = analytical * b0 * gamma_2pi;
simulation = simulation * b0 * gamma_2pi;
end

% Field in ppm
if (strcmp(unit, 'ppm'))
analytical = analytical * 1e6;
simulation = simulation * 1e6;
end



%% Plots to compare with analytical
close all
figure;
plot(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),dim_without_buffer(2)), ...
    squeeze(analytical(sectionx, sectiony, :)), ...
    'LineWidth',5,'Color',[0.6118    0.8824    1.0000],'LineStyle','-');
hold on
plot(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),dim_without_buffer(2)), ...
    squeeze(simulation(sectionx, sectiony, :)), ...
    'LineWidth',2,'Color',[0, 0.4470, 0.7410]);

hold off
xlabel('z position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
legend('Analytical', 'Simulation');
title(sprintf('Magnetic field %s [%s] in the %s phantom along z axis with susin=%0.2e and susout=%0.2e', field, unit, phantom, susin, susout))
grid on



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotting = 0;
if plotting == 1
    close all
    figure;
    plot(squeeze(dbz_analytical_ppm(sectionx, :, sectionz)));
    hold on
    plot(squeeze(dBz_map_ppm(sectionx, :, sectionz)));
    hold off
    xlabel('grid position')
    ylabel('dBz (ppm)')
    legend('Analytical', 'Simulation');
    title(sprintf('Field in the %s phantom along y axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))
    grid on
    
    figure;
    plot(squeeze(dbz_analytical_ppm(:, sectiony, sectionz)));
    hold on
    plot(squeeze(dBz_map_ppm(:, sectiony, sectionz)));
    hold off
    xlabel('grid position')
    ylabel('dBz (ppm)')
    legend('Analytical', 'Simulation');
    title(sprintf('Field in the %s phantom along x axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))
    grid on
    
    figure;
    plot(linspace(-dim_without_buffer(2)*res(2),dim_without_buffer(2)*res(2),dim_without_buffer(2)), ...
        squeeze(dbz_analytical_ppm(sectionx, sectiony, :)));
    hold on
    plot(linspace(-dim_without_buffer(2)*res(2),dim_without_buffer(2)*res(2),dim_without_buffer(2)), ...
        squeeze(dBz_map_ppm(sectionx, sectiony, :)));
    hold off
    xlabel('z position [mm]')
    ylabel('dBz [ppm]')
    legend('Analytical', 'simulation');
    title(sprintf('Field offset in the %s phantom along z axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))
    grid on


    %% Plot z-axis [ppm]
    close all
    figure;
    plot(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),dim_without_buffer(2)), ...
        squeeze(dbz_analytical_ppm(sectionx, sectiony, :)), ...
        'LineWidth',4,'Color',[0.6118    0.8824    1.0000],'LineStyle','-');
    hold on
    plot(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),dim_without_buffer(2)), ...
        squeeze(dBz_simulation_ppm(sectionx, sectiony, :)), ...
        'LineWidth',1,'Color',[0, 0.4470, 0.7410]);
    hold on
        plot(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),dim_without_buffer(2)), ...
        squeeze(Bz_analytical_ppm(sectionx, sectiony, :)), ...
        'LineWidth',4,'Color',	[1.0000    0.8471    0.4902],'LineStyle','-');
    hold on
    plot(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),dim_without_buffer(2)), ...
        squeeze(Bz_simulation_ppm(sectionx, sectiony, :)), ...
        'LineWidth',1,'Color',	[0.8500, 0.3250, 0.0980]);
    
    hold off
    xlabel('z position [mm]')
    ylabel('dBz [ppm]')
    legend('Analytical - offset', 'Simulation - offset','Analytical - demodulated', 'Simulation-demodulated');
    title(sprintf('Field offset [ppm] in the %s phantom along z axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))
    grid on

    %% Plot result at a section
    
    %   y section
    figure;
    subplot(1, 2, 1);
    plot_y_section = imagesc(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        squeeze(sus(:, sectiony, :))*1e6); colorbar; % colormap winter;
    title('Susceptibility distribution [ppm]');
    xlabel('z position [mm]')
    ylabel('x position [mm]')
    axis square

    subplot(1, 2, 2);
    imagesc(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        squeeze(dBz_simulation_ppm(:, sectiony, :))); colorbar;
    title('Simulation of the field variation [ppm]');
    xlabel('z position [mm]')
    ylabel('x position [mm]')
    axis square
    
    sgtitle(sprintf('section y = %d for %s with radius %0.2f mm', sectiony-dim_without_buffer(2)/2-1, 'sphere', radius))

    figure;
    subplot(1, 3, 1);
    imagesc(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        squeeze(sus(sectionx, :, :))*1e6); colorbar; % colormap winter;
    title('x');
    xlabel('z position [mm]')
    ylabel('y position [mm]')
    axis square
    
    subplot(1, 3, 2);
    imagesc(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        squeeze(sus(:, sectiony, :))*1e6); colorbar; % colormap winter;
    title('y');
    xlabel('z position [mm]')
    ylabel('x position [mm]')
    axis square

    subplot(1, 3, 3);
    imagesc(linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        linspace(-dim_without_buffer(2)/2*res(2),dim_without_buffer(2)/2*res(2),10), ...
        squeeze(sus(:, :,sectionz))*1e6); colorbar; % colormap winter;
    title('z');
    xlabel('y position [mm]')
    ylabel('x position [mm]')
    axis square
end    


    %% Plots to compare successive buffers
    % Run 3 simulations to store the results under the correct names. Make sure
    % the analytical matrix of the variation used has the maximum resolution.
    
    % figure;
    % plot(squeeze(anal(sectionx, sectiony, :)), 'LineWidth', 1);
    % hold on
    % plot(1:128, squeeze(dbz_128_z_cyly(sectionx, sectiony, :)), 'LineWidth', 1);
    % hold on
    % plot(1:128, squeeze( dbz_xz256_z_cyly(sectionx, sectiony, :)), 'LineWidth', 1);
    % hold on
    % plot(1:128, squeeze( dbz_xyz256_z_cyly(sectionx, sectiony, :)), 'LineWidth', 1);
    % hold off
    % xlabel('y position')
    % ylabel('dBz (ppm)')
    % legend('Analytical', 'simulation 128^3', 'simulation 256x128x256', 'simulation 256^3');
    % title(sprintf('Field in the %s phantom in ppm for theta = pi/2 with susin=1 and susout=0', phantom))
    % grid on

    %% Plot intermediate results
% figure; imagesc(squeeze((1/3-k_scaling_coeff(:, :, sectionz)))); colorbar;
% title('Multiplying volume in FT domain, x section')

% k_ifft = ifftshift(ifftn(ifftshift((1/3 - k_scaling_coeff))));
% 
% figure; imagesc(squeeze(k_ifft(sectionx, :, :))), colorbar;
% title('Ifft of the multiplying volume (spatial domain), x section')

