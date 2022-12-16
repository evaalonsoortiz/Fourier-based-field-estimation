%% Script to test subsample and the effect on resolution on the Zubal phantom
% Mathilde Dupouy made this script in 2022 for testing purposes. No further
% support is offered on the content of this script.

% run this script in the main folder (with FBFest.m) after adding the folder test_calc_dbz to the path
% add file 'zubal_modified.nii' to the main folder before running this script
% explanation on where to find this file is given in Zubal.m

% Many rests of precedent simulation are remaining but can inspire new
% experiments

%% 
%%    TEST SUB-SAMPLING
%%
%% Parameters
% dbz_path = 'zubal_downsamp/dbz_ZubalTest'
% b0 = 3; % [T]
%% generate susceptibility distribution for the modified Zubal phantom,
% previously downloaded
    % Properties of the phantom : 
    % dimensions 256x256x128
    % Resolution 1.1x1.1x1.4
zubal_sus_dist = Zubal('zubal_EAO.nii');
sus_nii = make_nii(zubal_sus_dist.volume);
save_nii(sus_nii, 'sus_zubal_EAO.nii')

dim_without_buffer = zubal_sus_dist.matrix;

sus = zubal_sus_dist.volume;

dim = 1 * dim_without_buffer;

%% Test imgaussfilt3
sigma = 0.5;
sus_filt_imgauss3 = imgaussfilt3(zubal_sus_dist.volume, sigma, 'FilterDomain', 'frequency', 'padding', 'symmetric');
% volume_gray = uint8(255*(mat2gray(sus_low_3D)));
% figure; montage(volume_gray); title('Test imgaussfilt3');

% Filtering in frequency domain using gaussmf
dim_pad = 2; % dim_pad pixels will be added in each direction
zsection =  dim_without_buffer(3) / 2 + 1;
xsection = dim_without_buffer(1) / 2 + 1;
sigma = [128, 128, 64];

[kx, ky, kz] = ndgrid(-dim_without_buffer(1)/2 - dim_pad + 1:dim_without_buffer(1)/2 + dim_pad, -dim_without_buffer(2)/2 - dim_pad + 1:dim_without_buffer(2)/2 + dim_pad, -dim_without_buffer(3)/2 - dim_pad + 1:dim_without_buffer(3)/2 + dim_pad);
gauss_3D = gaussmf(kx, [sigma(1), 0]) .* gaussmf(ky, [sigma(2), 0]) .* gaussmf(kz, [sigma(3), 0]);

figure; imagesc(squeeze(gauss_3D(:, :, zsection))); colorbar; title('3D gaussian visualisation');
volume_gray = uint8(255*mat2gray(gauss_3D));
figure; imagesc(squeeze(abs(ifftshift(ifftn(ifftshift(gauss_3D(:, :, 64))))))); colorbar; title('ifft of the 3D gaussian');

% Zubal FFT and filtering
zubal_fft = fftshift(fftn(fftshift(padarray(zubal_sus_dist.volume, [dim_pad dim_pad dim_pad], zubal_sus_dist.volume(1, 1, 1), 'symmetric'))));
zubal_filt_fft = zubal_fft .* gauss_3D;
zubal_filt = ifftshift(ifftn(ifftshift(zubal_filt_fft)));

zubal_filt = zubal_filt(dim_pad + 1: end - dim_pad, dim_pad + 1: end - dim_pad, dim_pad + 1: end - dim_pad);

figure; 
volume_gray_filt = uint8(255*mat2gray(real(zubal_filt)));
montage(volume_gray_filt); 
title('Filtered susceptibility')

figure; 
volume_gray_filt = uint8(255*mat2gray(real(sus)));
montage(volume_gray_filt); 
title('Initial susceptibility')

figure; 
volume_gray_filt = uint8(255*mat2gray(abs(real(zubal_filt) - sus)));
montage(volume_gray_filt); 
title('Difference between initial and filtered susceptibility')


%% Experiment 2 : Comparing calculation between the initial susceptibility and the filtered one

dBz_obj_filt = FBFest('Zubal', zubal_filt, zubal_sus_dist.image_res, dim_without_buffer, sus(1, 1, 1), dim);
dBz_obj = FBFest('Zubal', sus, zubal_sus_dist.image_res, dim_without_buffer, sus(1, 1, 1), dim);
dBz_map_ppm_filt = real(dBz_obj_filt.volume * 1e6);
dBz_map_ppm = real(dBz_obj.volume * 1e6);

%% Plot field shifts

figure; 
subplot(1, 3, 1)
imagesc(squeeze(dBz_map_ppm(xsection, :, :))); colorbar;
title('Field shift on the initial susceptibility')
subplot(1, 3, 2)
imagesc(squeeze(dBz_map_ppm_filt(xsection, :, :))); colorbar;
title('Field shift on the filtered susceptibility')
subplot(1, 3, 3)
imagesc(abs(squeeze(dBz_map_ppm(xsection, :, :))-squeeze(dBz_map_ppm_filt(xsection, :, :)))); colorbar;
title('Absolute difference')
sgtitle(sprintf('Field shifts (ppm) for the section at x=%u', xsection))

figure;
histogram(abs((dBz_map_ppm(:))-(dBz_map_ppm_filt(:))))
title('Absolute error between the field shifts on the initial and filtered susceptibility')
xlabel('Absolute error (ppm)')

%% Down-sampling and real part
zubal_downsamp = real(zubal_filt(1:2:end, 1:2:end, 1:2:end));

%% Filtering and down sampling using sub_sample_3D

zubal_filt = sub_sample_3D(zubal_sus_dist.volume, [2, 2, 2]);

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
