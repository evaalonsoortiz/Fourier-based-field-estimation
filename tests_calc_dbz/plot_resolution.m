
analytical = importdata("results/resolution_test/analytical_res01.mat");
simulation1 = importdata("results/resolution_test/analytical_res1.mat");
simulation05 = importdata("results/resolution_test/analytical_res2.mat");

%analytical = importdata("results/resolution_test/analytical_256_res05.mat");
%simulation1 = importdata("results/resolution_test/simulation_256_res1.mat");
%simulation05 = importdata("results/resolution_test/simulation_256_res05.mat");
%simulation2 = importdata("results/resolution_test/simulation_256_res2.mat");
res1 = 1;
res05 = 0.5;
res2 = 2;

close all
figure;
plot(linspace(-dim_without_buffer(2)/2*0.1,dim_without_buffer(2)/2*0.1,dim_without_buffer(2)), ...
    squeeze(analytical(sectionx, sectiony, :)), ...
    'LineWidth',3,'Color',[0.65 0.65 0.65],'LineStyle',':');

hold on
plot(linspace(-dim_without_buffer(2)/2*2,dim_without_buffer(2)/2*2,dim_without_buffer(2)), ...
    squeeze(simulation05(sectionx, sectiony, :)), ...
    'LineWidth',1,'Color','r');

hold on
plot(linspace(-dim_without_buffer(2)/2,dim_without_buffer(2)/2,dim_without_buffer(2)), ...
    squeeze(simulation1(sectionx, sectiony, :)), ...
    'LineWidth',1,'Color','b');

%hold on
%plot(linspace(-dim_without_buffer(2)/2*res2,dim_without_buffer(2)/2*res2,dim_without_buffer(2)), ...
%    squeeze(simulation2(sectionx, sectiony, :)), ...
%    'LineWidth',1,'Color',[0, 0.5, 0]);

hold off
xlabel('z position [mm]')
ylabel(sprintf('Field %s [%s]', field, unit))
%legend('Analytical', 'Simulation res = 0.5mm', 'Simulation res = 1mm', 'Simulation res = 2mm');
legend('Analytical res 0.1mm', 'Analytical res 2mm', 'Analytical res 1mm');
title(sprintf('Magnetic field %s [%s] in the %s phantom along z axis with susin=%0.2e and susout=%0.2e', field, unit, phantom, susin, susout))
xlim([-dim_without_buffer(2)/2*res05 dim_without_buffer(2)/2*res05])
grid on