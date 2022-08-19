%% Script to test subsample and the effect on resolution on the Zubal phantom
% Many rests of precedent simulation are remaining but can inspire new
% experiments

%% 
%%    TEST BUFFER
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
% Add a buffer
sus = padarray(sus, (dim - dim_without_buffer) / 2, sus(1, 1, 1));
padDim =  (dim - dim_without_buffer) / 2;


%% 
%%    TEST SUB-SAMPLING
%%
%% Parameters
% dbz_path = 'zubal_downsamp/dbz_ZubalTest';
% b0 = 3; % [T]
% %% generate susceptibility distribution for the modified Zubal phantom,
% % previously downloaded
%     % Properties of the phantom : 
%     % dimensions 256x256x128
%     % Resolution 1.1x1.1x1.4
% phantom = 'zubal'
% switch(phantom)
%     case 'zubal'
%         zubal_sus_dist = Zubal('zubal_EAO.nii');
%         sus_nii = make_nii(zubal_sus_dist.volume);
%         save_nii(sus_nii, 'sus_zubal_EAO.nii')
% 
%         dim_without_buffer = zubal_sus_dist.matrix;
%         dim = 1 * dim_without_buffer;
% 
%         sus = zubal_sus_dist.volume;
%     case 'rect'
%         dim_without_buffer = [8, 8, 8];
%         dim = 1 * dim_without_buffer;
%         res = [1,1,1];
%         sus = zeros(dim);
%         sus(2:5, 3:5, 4:7) = 1;
% end
% 
% % Add a buffer
% padDim =  (dim - dim_without_buffer) / 2;
% sus = padarray(sus, padDim, sus(1, 1, 1));

%% First filters
% filter = 'gauss'; param = 0.5;
% sus_low = sub_sample_section(zubal_sus_dist.volume, 'z', 65, filter, param);
% 
% figure; imagesc(sus_low(57:75, 102:111)); colorbar;
% title(sprintf('down-sampled %s filtered susceptibility map with a parameter %0.2f', filter, param));
% 
% sus_low_3D = sub_sample_3D(zubal_sus_dist.volume, filter, param);
% 
% figure; imagesc(sus_low_3D(57:75, 102:111, 33)); colorbar;
% title(sprintf('3D down-sampled %s filtered susceptibility map with a parameter %0.2f', filter, param));
% 
% figure;
% volume_gray = uint8(255*(1 - mat2gray(sus_low_3D)));
% montage(volume_gray);
% title(sprintf('3D down-sampled %s filtered susceptibility map with a parameter %0.2f', filter, param));
% 
% % plot difference
% figure;
% diff1 = abs(zubal_sus_dist.volume(1:2:end, 1:2:end, 1:2:end) - sus_low_3D);
% volume_gray = uint8(255*mat2gray(diff1));
% montage(volume_gray);
% title(sprintf('Absolute difference between 3D down-sampled susceptibility map and %s filtered with a parameter %0.2f', filter, param));

% Zubal FFT
% zubal_fft = fftshift(fftn(fftshift(zubal_sus_dist.volume)));
%
% figure; imagesc(log(1 + abs(zubal_fft(:, :, 65)))); colorbar;
% title(sprintf('3D susceptibility map FFT (with log transformation)'));

%% Test imgaussfilt3
sigma = 0.5;
sus_filt_imgauss3 = imgaussfilt3(zubal_sus_dist.volume, sigma, 'FilterDomain', 'frequency', 'padding', 'symmetric');
% volume_gray = uint8(255*(mat2gray(sus_low_3D)));
% figure; montage(volume_gray); title('Test imgaussfilt3');

% Filtering in frequency domain using gaussmf

dim_pad = 2; % dim_pad pixels will be added in each direction
%sigma = [128, 128, 64];
zsection =  dim(3) / 2 + 1;
xsection = dim(1) / 2 + 1;
sigma = [128, 128, 64];

[kx, ky, kz] = ndgrid(-dim(1)/2 - dim_pad + 1:dim(1)/2 + dim_pad, -dim(2)/2 - dim_pad + 1:dim(2)/2 + dim_pad, -dim(3)/2 - dim_pad + 1:dim(3)/2 + dim_pad);
gauss_3D = gaussmf(kx, [sigma(1), 0]) .* gaussmf(ky, [sigma(2), 0]) .* gaussmf(kz, [sigma(3), 0]);

figure; imagesc(squeeze(gauss_3D(:, :, zsection))); colorbar; title('3D gaussian visualisation');
volume_gray = uint8(255*mat2gray(gauss_3D));
figure; montage(volume_gray); title('3D gaussian visualisation bis')
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
title('Difference')

zubal_filt = sub_sample_3D(zubal_sus_dist.volume, 2);
%% Experiment 2 : Comparing calculation between the initial susceptibility and the filtered one

dBz_obj_filt = FBFest(zubal_filt, zubal_sus_dist.image_res, dim, b0);
dBz_obj = FBFest(sus, zubal_sus_dist.image_res, dim, b0);
dBz_map_ppm_filt = real(dBz_obj_filt.volume * 1e6);
dBz_map_ppm = real(dBz_obj.volume * 1e6);

%%
figure; 
volume_gray_filt = uint8(255*mat2gray(dBz_map_ppm_filt));
montage(volume_gray_filt); 
title('Field shift calculated on filtered susceptibility')

figure; 
volume_gray_filt = uint8(255*mat2gray(dBz_map_ppm));
montage(volume_gray_filt); 
title('Field shift calculated on initial susceptibility')

figure; 
volume_gray_filt = uint8(255*mat2gray(abs(dBz_map_ppm_filt - dBz_map_ppm)));
montage(volume_gray_filt); 
title('Difference between the field shift on the filtered susceptibilityand the initial one')

%%
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
imagesc((squeeze(dBz_map_ppm(xsection, :, :))-squeeze(dBz_map_ppm_filt(xsection, :, :))) ./ squeeze(dBz_map_ppm(xsection, :, :)))
colorbar;
relativeError = (abs(dBz_map_ppm(:)-(dBz_map_ppm_filt(:))) ./ abs(dBz_map_ppm(:)));
error = abs(dBz_map_ppm(:)-dBz_map_ppm_filt(:));
sqrt(sum(((dBz_map_ppm(:))-(dBz_map_ppm_filt(:))).^2) / prod(dim))

figure;
histogram(abs((dBz_map_ppm(:))-(dBz_map_ppm_filt(:))))

%% Down-sampling and real part
zubal_downsamp = real(zubal_filt(1:2:end, 1:2:end, 1:2:end));

%% Diffs

%% 
%%    TEST FIELD SHIFT ESTIMATION ON SUBSAMPLED
%%
%% New susceptibility and add buffer

% sus = zubal_downsamp;
% dim_without_buffer = size(sus);
% dim = 2 * dim_without_buffer; %+ [0, 0, 128];
% % Add a buffer
% sus = padarray(sus, (dim - dim_without_buffer) / 2, sus(1, 1, 1));
% padDim =  (dim - dim_without_buffer) / 2;
% 
% %% Variation calculation
% % tic
% dBz_obj_downsamp = FBFest(sus, zubal_sus_dist.image_res, dim, b0);
% % toc
% dBz_map_ppm_downsamp = real(dBz_obj.volume * 1e6); %TODO remove real ? Loss of the y-translation
% 
% %Truncate :
% dBz_obj_downsamp.matrix = dim_without_buffer;
% dBz_obj_downsamp.volume = dBz_obj_downsamp.volume( padDim(1) + 1:dim_without_buffer(1) + padDim(1), padDim(2) + 1:dim_without_buffer(2) + padDim(2), padDim(3) + 1:dim_without_buffer(3) + padDim(3));
% % Save NIFTI image of the field shift
% dBz_obj_downsamp.save([dbz_path '_downsamp']);
% 
% %% Down sampling the true dbz map
% dBz_obj.matrix = dBz_obj.matrix / 2;
% dBz_obj.volume = dBz_obj.volume(1:2:end, 1:2:end, 1:2:end); % Down sampling the calculated field shift
% 


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
