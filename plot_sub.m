ana = load('results/2511/ana_subsample.mat').analytical_x;
sim = load('results/2511/simulation_notsub.mat').simulation;
sim_sub = load('results/2511/simulation_subsample.mat').simulation;

sectionx_sub = round(dim_without_buffer(2)/sub_factor / 2) + 1;
sectiony_sub = round(dim_without_buffer(1)/sub_factor / 2) + 1;

% only along z-axis
figure;
plot(linspace(-dim_without_buffer(3)/2 * res_analytical(3),dim_without_buffer(3)/2 * res_analytical(3),dim_without_buffer(3)), ...
    squeeze(ana(sectionx, sectiony, :)), ...
    'LineWidth',3,'Color',[0.55 0.55 0.55],'LineStyle',':');
hold on
plot(linspace(-dim_without_buffer(3)/2*res(3),dim_without_buffer(3)/2*res(3),dim_without_buffer(3)), ...
    squeeze(sim(sectionx, sectiony, :)), ...
    'LineWidth',1.5,'Color','r');
hold on
plot(linspace(-dim_without_buffer(3)/2*res(3),dim_without_buffer(3)/2*res(3),dim_without_buffer(3)/sub_factor), ...
    squeeze(sim_sub(sectionx_subsample, sectiony_subsample, :)), ...
    'LineWidth',1.5,'Color','b');

hold off
xlabel('z position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
legend('Analytical', 'Simulation','Subsampled simulation');
title(sprintf('Magnetic field %s [%s] in %s (r=%0.1f mm) along z-axis with sus diff=%0.2f ppm', field, unit, phantom, radius, sus_diff*1e6))
grid on