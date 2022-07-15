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

