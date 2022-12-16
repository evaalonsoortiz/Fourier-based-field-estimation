function [] = plot_along_axes(simulation, unit)

sectionx = round(size(simulation,2) / 2) + 1;
sectiony = round(size(simulation,1) / 2) + 1;
sectionz = round(size(simulation,3) / 2) + 1;

max_lim = max(simulation, [],'all') + 0.1 * max(simulation, [],'all');
min_lim = min(simulation, [],'all') + 0.1 * min(simulation, [],'all');

close all

% plot along x, y and z axis of the simulated and analytical field
figure;
subplot(1, 3, 1);
plot(linspace(-size(simulation,2)/2,size(simulation,2)/2, size(simulation,2)), ...
    squeeze(simulation(sectiony, :, sectionz)), ...
    'LineWidth',2,'Color','k');

ylim([min_lim max_lim])
xlabel('x voxel')
ylabel(sprintf('Field [%s]', unit))
title('Field along x-axis')
grid on

subplot(1, 3, 2);
plot(linspace(-size(simulation,1)/2,size(simulation,1)/2, size(simulation,1)), ...
    squeeze(simulation(:, sectionx, sectionz)), ...
    'LineWidth',2,'Color','k');

ylim([min_lim max_lim])
xlabel('y voxel')
ylabel(sprintf('Field [%s]', unit))
title('Field along y-axis')
grid on


subplot(1, 3, 3);
plot(linspace(-size(simulation,3)/2,size(simulation,3)/2, size(simulation,3)), ...
    squeeze(simulation(sectiony, sectionx, :)), ...
    'LineWidth',2,'Color','k');

ylim([min_lim max_lim])
xlabel('z voxel')
ylabel(sprintf('Field [%s]', unit))
title('Field along z-axis')
grid on
