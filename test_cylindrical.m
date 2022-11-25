cyl_test = Cylindrical([100,100,100] , [1,1,1], 10, pi/2, [1 0]);
y_section = cyl_test.volume(50,:,:);
x_section = cyl_test.volume(:,50,:);
z_section = cyl_test.volume(:,:,50);

figure;
subplot(1,3,1)
imagesc(linspace(-50,50),linspace(-50,50),squeeze(x_section)); 
title('x section')
axis square

subplot(1,3,2)
imagesc(linspace(-50,50),linspace(-50,50),squeeze(y_section)); 
title('y section')
axis square

subplot(1,3,3)
imagesc(linspace(-50,50),linspace(-50,50),squeeze(z_section)); 
title('z section')
axis square
