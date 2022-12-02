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

dim_pad = dim(1)/2;  % dim_pad pixels will be added in each direction TODO adapted to sigma ?

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
vol_lowRes1 = vol_filt(phase(1):factor(1):end, phase(2):factor(2):end, phase(3):factor(3):end);
vol_lowRes2 = vol_filt(phase(1)+1:factor(1):end, phase(2)+1:factor(2):end, phase(3)+1:factor(3):end);
vol_lowRes = (vol_lowRes1 + vol_lowRes2) ./2;
end