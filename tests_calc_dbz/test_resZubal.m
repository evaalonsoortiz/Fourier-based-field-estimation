%% Script to test subsample and the effect on resolution on the Zubal phantom


% generate susceptibility distribution for the modified Zubal phantom
zubal_sus_dist = Zubal('zubal_EAO.nii');

% Properties of the phantom : dimensions 256x256x128
% Resolution 1.1x1.1x1.4
filter = 'gauss'; param = 0.5;
sus_low = sub_sample_section(zubal_sus_dist.volume, 'z', 65, filter, param);

figure; imagesc(sus_low(57:75, 102:111)); colorbar;
title(sprintf('down-sampled %s filtered susceptibility map with a parameter %0.2f', filter, param));

sus_low_3D = sub_sample_3D(zubal_sus_dist.volume, filter, param);

figure; imagesc(sus_low_3D(57:75, 102:111, 33)); colorbar;
title(sprintf('3D down-sampled %s filtered susceptibility map with a parameter %0.2f', filter, param));

figure;
volume_gray = uint8(255*(1 - mat2gray(sus_low_3D)));
montage(volume_gray);
title(sprintf('3D down-sampled %s filtered susceptibility map with a parameter %0.2f', filter, param));

figure;
diff = abs(zubal_sus_dist.volume(1:2:end, 1:2:end, 1:2:end) - sus_low_3D);
volume_gray = uint8(255*mat2gray(diff));
montage(volume_gray);
title(sprintf('Absolute difference between 3D down-sampled susceptibility map and %s filtered with a parameter %0.2f', filter, param));

%% Useful commands :

% Gray scaling / Inverse gray scaling
% volume_gray = ceil((zubal_sus_dist.volume - min(unique(zubal_sus_dist.volume)))*255/(max(unique(zubal_sus_dist.volume)) - min(unique(zubal_sus_dist.volume))));
% volume_gray = 255 - ceil((zubal_sus_dist.volume - min(unique(zubal_sus_dist.volume)))*255/(max(unique(zubal_sus_dist.volume)) - min(unique(zubal_sus_dist.volume))));
% %OR
% volume_gray = uint8(255*mat2gray(zubal_sus_dist.volume));
% volume_gray = uint8(255*(1 - mat2gray(zubal_sus_dist.volume)));

% To visualize all the slices of a 3D volume
% montage(volume_gray, 'DisplayRange', [0, 255])
% %OR
% montage(volume_gray)