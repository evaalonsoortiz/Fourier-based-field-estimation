x = 50;
test = zeros(x,x,x);
test(round(x/3),round(x/2):round(x/2)+10,:) = 100; % 

theta = 0;

t = [cos(theta)  0      sin(theta)   0
  0             1              0     0
  -sin(theta)    0       cos(theta)   0
  0             0              0     1];


tform = affine3d(t);

test_rotated = imwarp(test,tform);
%test_rotated = imrotate3(test, 90, [0 1 0]);


y_section = test_rotated(round(x/3),:,:);
x_section = test_rotated(:,round(x/2),:);
z_section = test_rotated(:,:,2*round(x/3));

%volshow(test_rotated)
close all

figure;
subplot(1,3,1)
imagesc(linspace(x/2,x/2),linspace(-x/2, x/2),squeeze(x_section)); 
title('x section')
axis square

subplot(1,3,2)
imagesc(linspace(x/2,x/2),linspace(-x/2, x/2),squeeze(y_section)); 
title('y section')
axis square

subplot(1,3,3)
imagesc(linspace(x/2,x/2),linspace(-x/2, x/2),squeeze(z_section)); 
title('z section')
axis square