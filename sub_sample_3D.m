function vol_lowRes = sub_sample_3D(volume, filter, param)
% Down sample the matrix by two with a previous low-pass filter.
%
% _SYNTAX_
% 
% vol_lowRes = sub_sample(matrix, filter, param)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    volume
%      A 3D matrix
%    filter
%       A string representing the low-pass filter wanted. It can be
%       'no-filter', 'average', 'gauss', 'median'
%    param
%       A parameter fot the filter : anything for 'no-filter', the kernel
%       side size for 'average' and 'median, sigma for 'gauss'.
%
% _OUTPUTS_
%
%   mat_lowRes
%     The down sampled matrix. Its size is two times lower than the matrix
%     size in each direction.
%

%% Low-pass filtering
switch(filter)
    case 'no-filter'
    case 'average'      
        volume = convn(volume, 1/param^2 * ones(param), 'same');
    case 'gauss'
        volume = imgaussfilt3(volume, param, 'padding', 'symmetric');
    case 'median'
        volume = medfilt3(volume, [param, param, param]);
    otherwise
end

%% Down sampling

vol_lowRes = volume(1:2:end, 1:2:end, 1:2:end);

end