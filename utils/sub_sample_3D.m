function vol_lowRes = sub_sample_3D(volume, factor, sigmaFilter)
% Down sample the matrix by two with a previous low-pass filter.
%
% _SYNTAX_
% 
% vol_lowRes = sub_sample_3D(phantomName, volume, factor, sigmaFilter)

% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    volume
%      A 3D matrix.
%    factor
%       An array of 3 elements representing the down-sampling factor in the three
%       directions (usually [2, 2, 2]).
%    sigmaFilter - optional
%       The cut-off frequency for the filter
%
% _OUTPUTS_
%
%   mat_lowRes
%     The down sampled matrix.
%

dim = size(volume);

%% Low-pass filtering
% Determination of sigma (equivalent to the cutoff frequency) if not passed
% in the parameters
if (nargin < 3)
    % To avoid aliasing, the cut-off frequency has to have fmax in the new
    % k-space view
    sigmaFilter = dim ./ factor;  
end

dim_pad = 0; %dim(1)/2;  % dim_pad pixels will be added in each direction TODO adapted to sigma ?

% Gaussian filter
%-dim(1)/2 - dim_pad  + 1:dim(1)/2 + dim_pad
[kx, ky, kz] = ndgrid(...
    -dim(1)/2 - dim_pad + 1:dim(1)/2 + dim_pad, ...
    -dim(2)/2 - dim_pad + 1:dim(2)/2 + dim_pad,...
    -dim(3)/2 - dim_pad + 1:dim(3)/2 + dim_pad);


gauss_3D = gaussmf(kx, [sigmaFilter(1), 0]) .* gaussmf(ky, [sigmaFilter(2), 0]) .* gaussmf(kz, [sigmaFilter(3), 0]);

% FFT with padding, filtering and IFFT
vol_fft = fftn(padarray(volume, [dim_pad dim_pad dim_pad], volume(1, 1, 1), 'symmetric'));
vol_filt = real(ifftn(vol_fft .* fftshift(gauss_3D)));
% Truncate padding
vol_filt = vol_filt(dim_pad + 1: end - dim_pad, dim_pad + 1: end - dim_pad, dim_pad + 1: end - dim_pad);

%% Down sampling
phase = rem(floor(dim / 2)+1, factor); % phase is chosen to be sure to keep the center at the center
phase(phase == 0) = factor(phase == 0);

% note: this subsampling is only symmetrical for factor 2, the interpolation scheme becomes a lot more difficult for larger factors
% because to stay symmetrical we need an array with dimension that follows mod(dimension, factor) = 1 
% (examples: for 2 this can be 3, 5, 7, ..., for 3 this can be 4, 7, 10, ..., for 4 this can be 5, 9, 13, ...)

if mod(factor(1),2) == 0
    % first interpolate in 3 dimensions if the factor is even
    vol_interp = interp3(vol_filt);

    % then take the interpolated values to create a matrix with odd dimensions
    vol_odd = vol_interp(2:2:end-1, 2:2:end-1, 2:2:end-1);

    % subsample
    vol_lowRes = vol_odd(1:factor(1):end, 1:factor(1):end, 1:factor(1):end);


else
    vol_lowRes = vol_filt(1:factor(1):end, 1:factor(1):end, 1:factor(1):end);
end
end