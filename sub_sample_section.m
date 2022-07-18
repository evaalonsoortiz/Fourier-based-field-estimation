function mat_lowRes = sub_sample_section(matrix, sectionName, sectionIndex, filter, param)
% Cut the matrix at the desired section and down sample the 2D section by
% two with a previous low-pass filter.
%
% _SYNTAX_
% 
% mat_lowRes = sub_sample(matrix, filter, param)
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    matrix
%      A 3D matrix
%    sectionName
%       The considered section in the 3D matrix. Either 'x', 'y' or 'z'.
%    sectionIndex
%       The index of the section to slice.
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
switch sectionName
    case 'x'
        mat_section = squeeze(matrix(sectionIndex, :, :));
    case 'y'
        mat_section = squeeze(matrix(:, sectionIndex, :));

    case 'z'
        mat_section = squeeze(matrix(:, :, sectionIndex));
    otherwise
end

%% Low-pass filtering
switch(filter)
    case 'no-filter'
    case 'average'      
        mat_section = conv2(mat_section, 1/param^2 * ones(param), 'same');
    case 'gauss'
        mat_section = imgaussfilt(mat_section, param, 'padding', 'symmetric');
    case 'median'
        mat_section = medfilt2(mat_section, [param, param]);
    otherwise
end

%% Down sampling

mat_lowRes = mat_section(1:2:end, 1:2:end);

end