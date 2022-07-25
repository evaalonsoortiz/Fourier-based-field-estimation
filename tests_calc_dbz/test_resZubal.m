%% Script to test subsample and the effect on resolution on the Zubal phantom

%% 
%%    TEST BUFFER
%%
%% Parameters
dbz_path = 'zubal_xyzbuffers/dbz_ZubalTest_buffer_trunc'
b0 = 3; % [T]
%% generate susceptibility distribution for the modified Zubal phantom,
% previously downloaded
    % Properties of the phantom : 
    % dimensions 256x256x128
    % Resolution 1.1x1.1x1.4
zubal_sus_dist = Zubal('zubal_EAO.nii');
sus_nii = make_nii(zubal_sus_dist.volume);
save_nii(sus_nii, 'sus_zubal_EAO.nii')

dim_without_buffer = zubal_sus_dist.matrix;

n = 40;
z_dims = 2 * ceil(linspace(0, 256, n)); % Be sure to have even and integers values
xy_dims = 2 * ceil(linspace(0, 128, n)); % Be sure to have even and integers values

% % Measures
% it_diffs = zeros(1, n); % store the iterative quadratique differences
% zsection = 65;
% xsection = 129;
% it_diffs_sec = zeros(256, 128, n-1);
% glob_diffs = zeros(1, n); %store the quadratique differences with the last 'better' calc
% glob_diffs_sec = zeros(256, 128, n);
% volumes_trunc = cell(1, n); % store the fields shift (ppm)
% temps = zeros(1, n);
% mean_value = zeros(1, n);
% 
% best = bestxy256_z512;
% 
% for zi = 1:n 
%     disp(z_dims(zi))
%     
%     sus = zubal_sus_dist.volume;
%     
%     dim = dim_without_buffer + [xy_dims(zi), xy_dims(zi), z_dims(zi)];
%     % Add a buffer
%     sus = padarray(sus, (dim - dim_without_buffer) / 2, sus(1, 1, 1), 'post');
%     sus = padarray(sus, (dim - dim_without_buffer) / 2, sus(1, 1, 1), 'pre');
%     padDim =  (dim - dim_without_buffer) / 2;
% 
% 
%     %% Variation calculation
%     % tic
%     t0 = cputime;
%     dBz_obj = FBFest(sus, zubal_sus_dist.image_res, dim, b0);
%     t1 = cputime;
% 
%     %Truncate :
%     dBz_obj.matrix = dim_without_buffer;
%     dBz_obj.volume = dBz_obj.volume( padDim(1) + 1:dim_without_buffer(1) + padDim(1), padDim(2) + 1:dim_without_buffer(2) + padDim(2), padDim(3) + 1:dim_without_buffer(3) + padDim(3));
%     % Save NIFTI image of the field shift
%     dBz_obj.save(sprintf('%s_y%u_z%u', dbz_path, xy_dims(zi), z_dims(zi)));
%     dBz_map_ppm = real(dBz_obj.volume * 1e6); %TODO remove real ? Loss of the y-translation
%     
%     volumes_trunc{zi} = dBz_map_ppm;
%     temps(zi) = t1-t0;
%     
%     if (zi > 1)
%         it_diffs(zi) = sum((dBz_map_ppm - volumes_trunc{zi - 1}) .^2, 'all')/prod(dim_without_buffer);
%         it_diffs_sec(:, :, zi-1) = squeeze((dBz_map_ppm(xsection, :, :) - volumes_trunc{zi - 1}(xsection, :, :)) .^2);
%     else
%         it_diffs(1) = sum((dBz_map_ppm - zeros(dim_without_buffer)) .^2, 'all')/prod(dim_without_buffer);
%     end
%     
%     glob_diffs(zi) =  sum((dBz_map_ppm - best) .^2, 'all')/prod(dim_without_buffer);
%     glob_diffs_sec(:, :, zi) = squeeze((dBz_map_ppm(xsection, :, :) - best(xsection, :, :)) .^2);
%     
%     mean_value(zi) = mean(dBz_map_ppm(:));
% 
% end
% %%
% figure; plot(z_dims, it_diffs, '.-');
% hold on 
% plot(z_dims, glob_diffs, '.-');
% hold off
% legend('Between successive iterations', 'Between the current calculation and the last (xy256 z 512)')
% xlabel('Pixels added in the z direction')
% ylabel('Quadratique error (ppm^2)')
% title('Mean of the quadratic errors while the dimension of the x, y and z-buffer increases')
% 
% figure;
% volume_gray = uint8(255*mat2gray(it_diffs_sec));
% montage(volume_gray, 'Size', [5, 8])
% 
% 
% figure;
% volume_gray = uint8(255*mat2gray(glob_diffs_sec));
% montage(volume_gray, 'Size', [5, 8])
% title('Evolution of the difference between the calculated fields for increasing x, y and z-buffers on a section at x=129 and the result for the last (bigger) buffer')
% 
% %%
% figure; 
% subplot(1, 3, 1)
% imagesc(squeeze(volumes_trunc{end-1}(xsection, :, :))); colorbar;
% title('Field shift without a buffer')
% subplot(1, 3, 2)
% imagesc(squeeze(volumes_trunc{end}(xsection, :, :))); colorbar;
% title('Field shift without a buffer xy256 and z512)')
% subplot(1, 3, 3)
% imagesc(abs(squeeze(volumes_trunc{end}(xsection, :, :))-squeeze(volumes_trunc{end-1}(xsection, :, :)))); colorbar;
% title('Absolute difference')
% sgtitle(sprintf('Field shifts for the section at x=%u', xsection))
% 
% %%
% figure;
% plot(z_dims, mean_value)
% 
% figure; 
% yyaxis left
% plot(z_dims, glob_diffs, '.-');
% yyaxis right
% plot(z_dims, mean_value, '.-')
% legend('quadratic error with the last volume iterations', 'Mean value in the volume')
% xlabel('Pixels added in the z direction')
% ylabel('Quadratique error (ppm^2)')

%%
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

% %% 
% %%    TEST SUB-SAMPLING
% %%
% %% First filters
% % filter = 'gauss'; param = 0.5;
% % sus_low = sub_sample_section(zubal_sus_dist.volume, 'z', 65, filter, param);
% % 
% % figure; imagesc(sus_low(57:75, 102:111)); colorbar;
% % title(sprintf('down-sampled %s filtered susceptibility map with a parameter %0.2f', filter, param));
% % 
% % sus_low_3D = sub_sample_3D(zubal_sus_dist.volume, filter, param);
% % 
% % figure; imagesc(sus_low_3D(57:75, 102:111, 33)); colorbar;
% % title(sprintf('3D down-sampled %s filtered susceptibility map with a parameter %0.2f', filter, param));
% % 
% % figure;
% % volume_gray = uint8(255*(1 - mat2gray(sus_low_3D)));
% % montage(volume_gray);
% % title(sprintf('3D down-sampled %s filtered susceptibility map with a parameter %0.2f', filter, param));
% % 
% % % plot difference
% % figure;
% % diff1 = abs(zubal_sus_dist.volume(1:2:end, 1:2:end, 1:2:end) - sus_low_3D);
% % volume_gray = uint8(255*mat2gray(diff1));
% % montage(volume_gray);
% % title(sprintf('Absolute difference between 3D down-sampled susceptibility map and %s filtered with a parameter %0.2f', filter, param));
% 
% %% Zubal FFT
% % zubal_fft = fftshift(fftn(fftshift(zubal_sus_dist.volume)));
% %
% % figure; imagesc(log(1 + abs(zubal_fft(:, :, 65)))); colorbar;
% % title(sprintf('3D susceptibility map FFT (with log transformation)'));
% 
% %% Test imgaussfilt3
% sigma = 0.5;
% sus_filt_imgauss3 = imgaussfilt3(zubal_sus_dist.volume, sigma, 'FilterDomain', 'frequency', 'padding', 'symmetric');
% % volume_gray = uint8(255*(mat2gray(sus_low_3D)));
% % figure; montage(volume_gray); title('Test imgaussfilt3');
% 
% %% Test filtering in frequency domain using gaussmf
% sigma = [128, 128, 64];
% % sigma = [64, 64, 32];
% 
% [kx, ky, kz] = ndgrid(-127:128, -127:128, -63:64);
% gauss_3D = gaussmf(kx, [sigma(1), 0]) .* gaussmf(ky, [sigma(2), 0]) .* gaussmf(kz, [sigma(3), 0]);
% % 
% % figure; imagesc(squeeze(gauss_3D(:, :, 65))); colorbar; title('3D gaussian visualisation');
% % volume_gray = uint8(255*mat2gray(gauss_3D));
% % figure; montage(volume_gray); title('3D gaussian visualisation bis')
% % figure; imagesc(squeeze(abs(ifftshift(ifftn(ifftshift(gauss_3D(:, :, 64))))))); colorbar; title('ifft of the 3D gaussian');
% 
% % Zubal FFT and filtering
% zubal_fft = fftshift(fftn(fftshift(zubal_sus_dist.volume)));
% zubal_filt_fft = zubal_fft .* gauss_3D;
% zubal_filt = ifftshift(ifftn(ifftshift(zubal_filt_fft)));
% 
% % figure; 
% % volume_gray_filt = uint8(255*mat2gray(abs(zubal_filt)));
% % montage(volume_gray_filt); 
% 
% %% Down-sampling and real part
% zubal_downsamp = real(zubal_filt(1:2:end, 1:2:end, 1:2:end));
% 
% %% Diffs
% % figure;
% % diff2 = real(zubal_filt(1:2:end, 1:2:end, 1:2:end)) - zubal_sus_dist.volume(1:2:end, 1:2:end, 1:2:end);
% % volume_gray = uint8(255*mat2gray(abs(diff2)));
% % montage(volume_gray);
% % title(sprintf('Absolute difference between 3D down-sampled susceptibility map and Gaussian filtered with a STD 64x64x32'));
% 
% %% 
% %%    TEST FIELD SHIFT ESTIMATION ON SUBSAMPLED
% %%
% %% New susceptibility and add buffer
% 
% sus = zubal_downsamp;
% dim_without_buffer = size(sus);
% dim = 2 * dim_without_buffer; %+ [0, 0, 128];
% % Add a buffer
% sus = padarray(sus, (dim - dim_without_buffer) / 2, sus(1, 1, 1), 'post');
% sus = padarray(sus, (dim - dim_without_buffer) / 2, sus(1, 1, 1), 'pre');
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
% %% Comparison
% diff3 = abs(dBz_obj.volume - dBz_obj_downsamp.volume);
% 
% volume_gray = uint8(255*mat2gray(diff3));
% montage(volume_gray);
% title(sprintf('Absolute difference between 3D down-sampled field map (ppm) and Gaussian filtered with a STD 64x64x32'));
% 
% 
% figure; imagesc(squeeze(real(dBz_obj.volume(:, :, 33)))); colorbar;
% title('Section of the  3D down-sampled field map (ppm)')
% 
% figure; imagesc(squeeze(real(dBz_obj_downsamp.volume(:, :, 33)))); colorbar;
% title('Section of Gaussian filtered with a STD 64x64x32 and down-sampled field map (ppm)')
% 
% figure; imagesc(squeeze(diff3(:, :, 33))); colorbar;
% title('Section of the absolute difference between 3D down-sampled field map (ppm) and Gaussian filtered with a STD 64x64x32')
% store_disp{1} = squeeze(real(dBz_obj.volume(:, :, 33)));
% store_disp{2} = squeeze(real(dBz_obj_downsamp.volume(:, :, 33)));
% store_disp{3} = squeeze(diff3(:, :, 33));
% a = +imutils.displayExperiment(3, store_disp, [0], ...
%     {'3D down-sampled field map (ppm)', 'Gaussian filtered with a STD 40x40x20 and down-sampled field map (ppm)', 'Difference'}, '', 'Section of the  3D down-sampled field map (ppm)');
% 
% % max de la difference : max(diff3, [], 'all') is 1.0001e-05
% % mean de la difference : max(diff3, [], 'all') is 1.0921e-07
% % mean of true ppm : mean(abs(dBz_obj.volume),'all') is 2.3895e-06
% % mean of true ppm : mean(abs(dBz_obj_downsamp.volume),'all') is  2.3355e-06
% %% Comparison with different sigma
% sectionz = 50;
% store_disp{1} = squeeze(real(dBz_obj.volume(:, :, sectionz)));
% diff3 = abs(dBz_obj.volume - dBz_obj.volume);
% store_disp{5} = squeeze(diff3(:, :, sectionz));
% 
% store_disp{2} = squeeze(real(vol_ds_323216(:, :, sectionz)));
% diff3 = abs(dBz_obj.volume - vol_ds_323216);
% store_disp{6} = squeeze(diff3(:, :, sectionz));
% 
% store_disp{3} = squeeze(real(vol_ds_646432(:, :, sectionz)));
% diff3 = abs(dBz_obj.volume - vol_ds_646432);
% store_disp{7} = squeeze(diff3(:, :, sectionz));
% 
% store_disp{4} = squeeze(real(vol_ds_12812864(:, :, sectionz)));
% diff3 = abs(dBz_obj.volume - vol_ds_12812864);
% store_disp{8} = squeeze(diff3(:, :, sectionz));
% 
% a = +imutils.displayExperiment(2, store_disp, [0, 0.5, 1, 2], ...
%     {'Gaussian filtered down-sampled field map (ppm), STD (64x64x32)x', 'Difference'}, 'factor', 'Section of the  3D down-sampled field map (ppm)')
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
