%% Script to test subsample and the effect on resolution on the Zubal phantom
%% Parameters
dbz_path = 'dbz_ZubalTest_bufferxy512z256_trunc'

%% generate susceptibility distribution for the modified Zubal phantom,
% previously downloaded
    % Properties of the phantom : 
    % dimensions 256x256x128
    % Resolution 1.1x1.1x1.4
zubal_sus_dist = Zubal('zubal_EAO.nii');
b0 = 3; % [T]

sus = zubal_sus_dist.volume;
dim_without_buffer = zubal_sus_dist.matrix;
dim = 2 * dim_without_buffer; %+ [0, 0, 128];
% Add a buffer
sus = padarray(sus, (dim - dim_without_buffer) / 2, sus(1, 1, 1), 'post');
sus = padarray(sus, (dim - dim_without_buffer) / 2, sus(1, 1, 1), 'pre');
padDim =  (dim - dim_without_buffer) / 2;


%% Variation calculation
% tic
dBz_obj = FBFest(sus, zubal_sus_dist.image_res, dim, b0);
% toc
dBz_map_ppm = real(dBz_obj.volume * 1e6); %TODO remove real ? Loss of the y-translation

%Truncate :
dBz_obj.matrix = dim_without_buffer;
dBz_obj.volume = dBz_obj.volume( padDim(1) + 1:dim_without_buffer(1) + padDim(1), padDim(2) + 1:dim_without_buffer(2) + padDim(2), padDim(3) + 1:dim_without_buffer(3) + padDim(3));
% Save NIFTI image of the field shift
dBz_obj.save(dbz_path);

% tic
%     %%---------------------------------------------------------------------- %%
%     %% Define constants
%     %%---------------------------------------------------------------------- %%
% 
%     % k-space window
%     k_max = 1./(2.*res);
%     interval = 2 * k_max ./ dim;
% 
%     % define k-space grid
%     [kx,ky,kz] = ndgrid(-k_max(1):interval(1):k_max(1) - interval(1),-k_max(2):interval(2):k_max(2) - interval(2),-k_max(3):interval(3):k_max(3) - interval(3));
% 
%     %%---------------------------------------------------------------------- %%
%     %% Compute Bdz
%     %%---------------------------------------------------------------------- %%
% 
%     % compute the fourier transform of the susceptibility distribution
%     FT_chi = fftshift(fftn(fftshift(sus)));
% 
%     % calculate the scaling coefficient 'kz^2/k^2'
%     k2 = kx.^2 + ky.^2 + kz.^2;
%     k_scaling_coeff = kz.^2./k2;
%     k_scaling_coeff(k2 == 0) = 0;
% 
%     % compute Bdz (the z-component of the magnetic field due to a
%     % sphere, relative to B0) FFT. The region of interest is
%     % assumed to be surrounded by a region of infinite extent whose
%     % susceptibility is equal to the susceptibility on the origin
%     % corner of the matrix.
%     bdzFFT = b0 * (1/3 - k_scaling_coeff).*FT_chi;
%     bdzFFT(k2 == 0) = b0 * sus(1, 1, 1) * prod(dim) / 3;
%     dbz_volume = ifftshift(ifftn(ifftshift(bdzFFT)));
% toc

% Conversion into ppm
%dBz_map_ppm = real(volume * 1e6); %TODO remove real ? Loss of the y-translation
%dBz_map_ppm = dbz_volume * 1e6; %TODO remove real ? Loss of the y-translation

% Properties of the phantom : dimensions 256x256x128
% Resolution 1.1x1.1x1.4

%% First filters
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

% plot difference
figure;
diff1 = abs(zubal_sus_dist.volume(1:2:end, 1:2:end, 1:2:end) - sus_low_3D);
volume_gray = uint8(255*mat2gray(diff1));
montage(volume_gray);
title(sprintf('Absolute difference between 3D down-sampled susceptibility map and %s filtered with a parameter %0.2f', filter, param));

%% Zubal FFT
% zubal_fft = fftshift(fftn(fftshift(zubal_sus_dist.volume)));
%
% figure; imagesc(log(1 + abs(zubal_fft(:, :, 65)))); colorbar;
% title(sprintf('3D susceptibility map FFT (with log transformation)'));

%% Test imgaussfilt3
sigma = 0.5;
sus_filt_imgauss3 = imgaussfilt3(zubal_sus_dist.volume, sigma, 'FilterDomain', 'frequency', 'padding', 'symmetric');
% volume_gray = uint8(255*(mat2gray(sus_low_3D)));
% figure; montage(volume_gray); title('Test imgaussfilt3');

%% Test filtering in frequency domain using gaussmf
sigma = [64, 64, 32];

[kx, ky, kz] = ndgrid(-127:128, -127:128, -63:64);
gauss_3D = gaussmf(kx, [sigma(1), 0]) .* gaussmf(ky, [sigma(2), 0]) .* gaussmf(kz, [sigma(3), 0]);

figure; imagesc(squeeze(gauss_3D(:, :, 65))); colorbar; title('3D gaussian visualisation');
volume_gray = uint8(255*mat2gray(gauss_3D));
figure; montage(volume_gray); title('3D gaussian visualisation bis')
figure; imagesc(squeeze(abs(ifftshift(ifftn(ifftshift(gauss_3D(:, :, 64))))))); colorbar; title('ifft of the 3D gaussian');

% Zubal FFT and filtering
zubal_fft = fftshift(fftn(fftshift(zubal_sus_dist.volume)));
zubal_filt_fft = zubal_fft .* gauss_3D;
zubal_filt = ifftshift(ifftn(ifftshift(zubal_filt_fft)));

figure; 
volume_gray_filt = uint8(255*mat2gray(abs(zubal_filt)));
montage(volume_gray_filt); 

%% Down-sampling and real part
zubal_downsamp = real(zubal_filt(1:2:end, 1:2:end, 1:2:end));

%% Diffs
figure;
diff2 = real(zubal_filt(1:2:end, 1:2:end, 1:2:end)) - zubal_sus_dist.volume(1:2:end, 1:2:end, 1:2:end);
volume_gray = uint8(255*mat2gray(abs(diff2)));
montage(volume_gray);
title(sprintf('Absolute difference between 3D down-sampled susceptibility map and Gaussian filtered with a STD 64x64x32'));

%% Test of field shift map determination



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
