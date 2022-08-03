%% Script to test adding buffers on the Zubal phantom

%% 
%%    TEST BUFFER
%%
%% Parameters
dbz_path = 'zubal_downsamp/dbz_ZubalTest'
b0 = 3; % [T]
%% generate susceptibility distribution for the modified Zubal phantom,
% previously downloaded
    % Properties of the phantom : 
    % dimensions 256x256x128
    % Resolution 1.1x1.1x1.4
    
zubal_sus_dist = Zubal('zubal_EAO.nii');
% sus_nii = make_nii(zubal_sus_dist.volume);
% save_nii(sus_nii, 'sus_zubal_EAO.nii')

dim_without_buffer = zubal_sus_dist.matrix;
sus = zubal_sus_dist.volume;
sus(sus == sus(1, 1, 1)) = sus(1,1,1);
sus_ext = sus(1,1,1);
sus_pad = 0;
sus = sus - sus_ext;

fprintf('Check if the susceptibility %0.4e is the external susceptibility.\n', sus_ext)
% dim = 1 * dim_without_buffer;
% Add a buffer
% padDim =  (dim - dim_without_buffer) / 2;
% sus = padarray(sus, padDim, sus_ext);

zubal_sus_dist.volume = sus;

%% Variation calculation "template"
% % tic
% disp(size(sus))
% t0 = cputime;
% dBz_obj = FBFest(sus, zubal_sus_dist.image_res, dim, b0);
% t1 = cputime;
% 
% %Truncate :
% dBz_obj.matrix = dim_without_buffer;
% dBz_obj.volume = dBz_obj.volume( padDim(1) + 1:dim_without_buffer(1) + padDim(1), padDim(2) + 1:dim_without_buffer(2) + padDim(2), padDim(3) + 1:dim_without_buffer(3) + padDim(3));
% % Save NIFTI image of the field shift
% %dBz_obj.save(sprintf('%s_y%u_z%u', dbz_path, xy_dims(zi), z_dims(zi)));
% dBz_map_ppm = real(dBz_obj.volume * 1e6); %TODO remove real ? Loss of the y-translation
    
%% Experiment 1 : Variation calculation with different buffers

n = 1;
z_dims = (2 * ceil(linspace(0, 0, n))); % Be sure to have even and integer values
xy_dims = (2 * ceil(linspace(128, 128, n))); % Be sure to have even and integer values

%% Measures
it_diffs = zeros(1, n); % store the iterative quadratique differences
glob_diffs = zeros(1, n); %store the quadratique differences with the last 'better' calc
zsection = 65;
xsection = 129;
it_diffs_sec = zeros(256, 128, n);
glob_diffs_sec = zeros(256, 128, n);
volumes_trunc = cell(1, n); % store the fields shift (ppm)
volumes_not_trunc = cell(1, n);
glob_diffs_sec_notTrunc = cell(1, n);
mean_value = zeros(1, n);

time = zeros(1, n);

best = zeros(dim_without_buffer);
%best = best_y512_air_Corr; % HERE
%best = best_z512; % HERE
%best = best_z512_minusSusAir_susext0;

for zi = 1:n 
    disp(z_dims(zi))
    
    sus = zubal_sus_dist.volume;
    
     % dim = dim_without_buffer + [xy_dims(zi), xy_dims(zi), z_dims(zi)];
     % dim = dim_without_buffer + [xy_dims(zi), 0, 0];
     dim = dim_without_buffer + [0, 0, z_dims(zi)];

    % Add a buffer
    padDim =  (dim - dim_without_buffer) / 2;
    sus = padarray(sus, padDim, sus_pad);
    
    %% Variation calculation
    % tic
    t0 = cputime;
    dBz_obj = FBFest(sus, zubal_sus_dist.image_res, dim, sus_ext, b0);
    t1 = cputime;
    
    %%
    FT_chi = fftshift(fftn(fftshift(sus + sus_ext)));
    FT_s = fftshift(fftn(fftshift(sus)));
     % k-space window
    k_max = 1./(2.*zubal_sus_dist.image_res);
    interval = 2 * k_max ./ dim;
    % define k-space grid
    [kx,ky,kz] = ndgrid(-k_max(1):interval(1):k_max(1) - interval(1),-k_max(2):interval(2):k_max(2) - interval(2),-k_max(3):interval(3):k_max(3) - interval(3));
    % calculate the scaling coefficient 'kz^2/k^2'
    k2 = kx.^2 + ky.^2 + kz.^2;
    kernel = b0 * (1/3 - kz.^2 ./ k2);
    kernel(k2 == 0) = b0 / 3;
    bdzFFT = kernel .* FT_chi;
    
    N = prod(dim);
    Ni = sum(sus ~= 0, 'all');
    sum_sn = sum(sus, 'all');
    moy_val = N*b0 / 3 * (sus_ext + sum_sn / N);
    
    bdz = ifftshift(ifftn(ifftshift(bdzFFT)));
    
%%     
%     if (zi == 3 || zi==4 || zi==2 || zi==n)
%         volumes_not_trunc{zi} = real(dBz_obj.volume * 1e6);
%         glob_diffs_sec_notTrunc{zi} =  squeeze(real(dBz_obj.volume(xsection, :, :) * 1e6));
%         glob_diffs_sec_notTrunc{zi}(padDim(2) + 1:dim_without_buffer(2) + padDim(2), padDim(3) + 1:dim_without_buffer(3) + padDim(3)) = (glob_diffs_sec_notTrunc{zi}(padDim(2) + 1:dim_without_buffer(2) + padDim(2), padDim(3) + 1:dim_without_buffer(3) + padDim(3)) - squeeze(best(xsection, :, :))) .^2;
%     end
    %Truncate :
    dBz_obj.matrix = dim_without_buffer;
    dBz_obj.volume = dBz_obj.volume( padDim(1) + 1:dim_without_buffer(1) + padDim(1), padDim(2) + 1:dim_without_buffer(2) + padDim(2), padDim(3) + 1:dim_without_buffer(3) + padDim(3));
    % Save NIFTI image of the field shift
    dBz_obj.save(sprintf('%s_y%u_z%u', dbz_path, xy_dims(zi), z_dims(zi)));
    dBz_map_ppm = real(dBz_obj.volume * 1e6); %TODO remove real ? Loss of the y-translation
    
    volumes_trunc{zi} = dBz_map_ppm;
    time(zi) = t1-t0;
    
    if (zi > 1)
        it_diffs(zi) = sum((dBz_map_ppm - volumes_trunc{zi - 1}).^2, 'all') / prod(dim_without_buffer);
        it_diffs_sec(:, :, zi) = squeeze((dBz_map_ppm(xsection, :, :) - volumes_trunc{zi - 1}(xsection, :, :)).^2);
    end
    
    glob_diffs(zi) =  sum((dBz_map_ppm - best) .^2, 'all') / prod(dim_without_buffer);
    glob_diffs_sec(:, :, zi) = squeeze((dBz_map_ppm(xsection, :, :) - best(xsection, :, :)) .^2);
    
    mean_value(zi) = mean(dBz_map_ppm(:));

end
%%
yyaxis left
figure; plot(z_dims(2:end), it_diffs(2:end), '.-'); %HERE
hold on 
plot(z_dims, glob_diffs, '.-'); %HERE
ylabel('Quadratic error (ppm^2)')
hold on
yyaxis right
plot(z_dims, time);
ylabel('time (ms)')
hold off
legend('Between successive iterations', 'Between the current calculation and the last (z512)', 'execution time (ms)') %HERE
xlabel('Pixels added in the z direction') %HERE
grid on
grid minor
title('Mean of the quadratic errors while the dimension of the z-buffer increases and x and y buffers are large') %HERE
%title('Mean of the quadratic errors while the dimension of the x, y and z-buffer increases')



% figure;
% volume_gray = uint8(255*mat2gray(it_diffs_sec));
% montage(volume_gray, 'Size', [5, 8])

figure;
yyaxis left
plot(z_dims, glob_diffs, '.-');
ylabel('Quadratic difference (ppm^2)') %HERE
yyaxis right
plot(z_dims, mean_value, '.-');
ylabel('Mean value in the volume (ppm)')
xlabel('Pixels added in the x direction') %HERE
legend('Quadratic difference between current and last iteration (xy256 z512)', 'Mean value in the volume')
title('Quadratic difference and mean value while the z-buffer increases and x and y buffers are large') %HERE
grid on
grid minor

%%
figure;
imagesc(squeeze(best_z512_susextAir(xsection, :, :) - best_z512_minusSusAir_susextAir(xsection, :, :))); colorbar;

% 
% figure;
% imagesc(squeeze(abs(TF_sus_z512susextAir_0(xsection, :, :)))); colorbar;
% figure;
% imagesc(squeeze(abs(TF_sus_z512susextAir_minus_0(xsection, :, :)))); colorbar;
% figure;
% imagesc(squeeze(abs(TF_sus_z512susextAir_minus_0(xsection, :, :))- abs(TF_sus_z512susextAir_0(xsection, :, :)))); colorbar;
% 
% figure;
% imagesc(squeeze(abs(TF_sus_z512susextAir(xsection, :, :)))); colorbar;
% figure;
% imagesc(squeeze(abs(TF_sus_z512susextAir_minus(xsection, :, :)))); colorbar;
% figure;
% imagesc(squeeze(abs(TF_sus_z512susextAir_minus(xsection, :, :))- abs(TF_sus_z512susextAir(xsection, :, :)))); colorbar;
% 
% figure;
% plot(squeeze(abs(TF_sus_z512susextAir_minus(xsection, 129, :))- abs(TF_sus_z512susextAir(xsection, 129, :))))

%% 

% figure; 
% subplot(2, 4, 1);
% imagesc(squeeze(volumes_not_trunc{2}(xsection, :, :))); colorbar;
% subplot(2, 4, 2);
% imagesc(squeeze(volumes_not_trunc{3}(xsection, :, :))); colorbar;
% subplot(2, 4, 3);
% imagesc(squeeze(volumes_not_trunc{4}(xsection, :, :))); colorbar;
% subplot(2, 4, 4);
% imagesc(squeeze(volumes_not_trunc{end}(xsection, :, :))); colorbar;
% 
% subplot(2, 4, 5);
% imagesc(squeeze(glob_diffs_sec_notTrunc{2})); colorbar;
% subplot(2, 4, 6);
% imagesc(squeeze(glob_diffs_sec_notTrunc{3})); colorbar;
% subplot(2, 4, 7);
% imagesc(squeeze(glob_diffs_sec_notTrunc{4})); colorbar;
% subplot(2, 4, 8);
% imagesc(squeeze(glob_diffs_sec_notTrunc{end})); colorbar;

% %%
% 
% s3 = size(volumes_not_trunc{3});
% repet3 = repmat(volumes_not_trunc{3}, 3, 3, 3);
% %%
% figure;
% imagesc(squeeze(repet3(s3(1) + xsection, :, :))); colorbar;
% title('Field shift (ppm) with a z-buffer of 28 voxels at a xsection ')
% %%
% 
% sn = size(volumes_not_trunc{end});
% repetn = repmat(volumes_not_trunc{end}, 3, 3, 3);
% %%
% figure;
% imagesc(squeeze(repetn(sn(1) + xsection, :, :))); colorbar;
% title('Field shift (ppm) with a z-buffer of 512 voxels at a xsection ')

%% Calculation "outside" FBFest
% tic
%     %%---------------------------------------------------------------------- %%
%     %% Define constants
%     %%---------------------------------------------------------------------- %%
% 
%     % k-space window
%     k_max = 1./(2.*[1.1, 1.1, 1.3]);
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
% 
% % Conversion into ppm
% dBz_map_ppm = real(volume * 1e6); %TODO remove real ? Loss of the y-translation
% dBz_map_ppm = dbz_volume * 1e6; %TODO remove real ? Loss of the y-translation
% 
% % Properties of the phantom : dimensions 256x256x128
% % Resolution 1.1x1.1x1.4