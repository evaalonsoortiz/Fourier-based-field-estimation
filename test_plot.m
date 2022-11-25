dim_without_buffer = [256 256 256];
res_analytical = [1 1 1];
res = [1 1 1];

sus = Cylindrical([256 256 256], [1 1 1], 10, pi/2, [3e-6 -3e-6]).volume;
simulation = FBFest('cylinder',sus, [1 1 1], [256 256 256],sus(1,1,1), [256 256 256]).volume;

sectionx = round(dim_without_buffer(2) / 2) + 1;
sectiony = round(dim_without_buffer(1) / 2) + 1;
sectionz = round(dim_without_buffer(3) / 2) + 1;

figure;

plot(linspace(-dim_without_buffer(3)/2*res(3),dim_without_buffer(3)/2*res(3),dim_without_buffer(3)), ...
    squeeze(simulation(sectionx, sectiony, :)), ...
    'LineWidth',1.5,'Color','r');

hold off
xlabel('z position [mm]')
%ylabel(sprintf('Field %s [%s]', field, unit))
%legend('Analytical', 'Simulation');
%title(sprintf('Magnetic field %s [%s] in %s (r=%0.1f mm) along z-axis with sus diff=%0.2f ppm', field, unit, phantom, radius, sus_diff*1e6))
grid on