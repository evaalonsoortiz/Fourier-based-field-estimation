sphere128 = load('cyl_r48_dim128.mat');
sphere_an128 = load('cyl_analytical_r48_dim128.mat');
sphere256 = load('cyl_r48_dim256.mat');
sphere384 = load('cyl_r48_dim384.mat');
%sphere512 = load('sphere_r48_512.mat');
%spherebuffer128 = load('sphere128buffered2.mat');

sectiony_128 = round(size(sphere128.dBz_map_ppm,2) / 2) + 1;
sectionx_128 = round(size(sphere128.dBz_map_ppm,1) / 2) + 1;
sectiony_256 = round(size(sphere256.dBz_map_ppm,2) / 2) + 1;
sectionx_256 = round(size(sphere256.dBz_map_ppm,1) / 2) + 1;
sectiony_384 = round(size(sphere384.dBz_map_ppm,2) / 2) + 1;
sectionx_384 = round(size(sphere384.dBz_map_ppm,1) / 2) + 1;

%% Plots to compare with analytical
close all % close previous plots


figure;
x1 = -63:64;
x2 = -127:128;
x3 = -191:192;
x_128 = linspace(-1,1,numel(sphere128.dBz_map_ppm(sectionx_128, sectiony_128, :)));
x_256 = linspace(-1,1,numel(sphere256.dBz_map_ppm(sectionx_256, sectiony_256, :)));

plot(x,squeeze(sphere_an128.dbz_analytical_ppm(sectionx_128, sectiony_128, :)));
hold on
plot(x,squeeze(sphere128.dBz_map_ppm(sectionx_128, sectiony_128, :)));
hold on
plot(x2,squeeze(sphere256.dBz_map_ppm(sectionx_256, sectiony_256, :)));
hold on
plot(x3,squeeze(sphere384.dBz_map_ppm(sectionx_384, sectiony_384, :)));
hold off
xlabel('grid position')
ylabel('dBz (ppm)')
legend('Analytical', 'Simulation 128^3', 'Simulation 256^3', 'Simulation 384^3');
title(sprintf('Field in the %s phantom in ppm along z axis with susin=%0.2e and susout=%0.2e', phantom, susin, susout))
grid on


