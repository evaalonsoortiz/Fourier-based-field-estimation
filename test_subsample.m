cyl_test = Cylindrical([100,100,100] , [1,1,1], 20, pi/2, [1 0]);
cyl_volume = cyl_test.volume();
factor = 5;
cyl_volume_sub = sub_sample_3D(cyl_volume,factor*[1,1,1]);

y_section = cyl_volume(50,:,:);
x_section = cyl_volume(:,50,:);
z_section = cyl_volume(:,:,50);

y_section_sub = cyl_volume_sub(round(50/factor),:,:);
x_section_sub = cyl_volume_sub(:,round(50/factor),:);
z_section_sub = cyl_volume_sub(:,:,round(50/factor));

figure;
subplot(2,3,1)
imagesc(linspace(-50,50),linspace(-50,50),squeeze(x_section)); 
title('x section')
axis square

subplot(2,3,2)
imagesc(linspace(-50,50),linspace(-50,50),squeeze(y_section)); 
title('y section')
axis square

subplot(2,3,3)
imagesc(linspace(-50,50),linspace(-50,50),squeeze(z_section)); 
title('z section')
axis square

subplot(2,3,4)
imagesc(linspace(-50,50),linspace(-50,50),squeeze(x_section_sub)); 
title('x section')
axis square

subplot(2,3,5)
imagesc(linspace(-50,50),linspace(-50,50),squeeze(y_section_sub)); 
title('y section')
axis square

subplot(2,3,6)
imagesc(linspace(-50,50),linspace(-50,50),squeeze(z_section_sub)); 
title('z section')
axis square