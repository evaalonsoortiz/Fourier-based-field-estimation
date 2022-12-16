function [distROI, d8_matrix] = calc_dist_ROI(matrix)
% Mathilde Dupouy made this script in 2022 for buffer dimension calculations. No further
% support is offered on the content of this script.

% This function calculates the minimum 8-distance between the region of
% interest (ROI, defined as what is not the background) and the edge of the
% 3D matrix. The method is based on two spanning of the matrixes and at
% each point the distance at the current point is compared with the
% distances + 1 in the neighbours, where the neighbours are up and left
% then down and right, with the spanning of the matrix respectively from
% the upper left and the lower right.
%
% This function can be slow, there are certainly improvements to make.
%
% _SYNTAX_
% 
% [distROI, distMat] = calc_dist_ROI(matrix);
%
% _DESCRIPTION_
%
% _INPUT ARGUMENTS_
%
%    matrix
%      A 3D matrix with an element surrounded by a uniform background. In
%      particular, the value in the (1, 1, 1) corner has to be part of the
%      background.
%
% _OUTPUTS_
%
%   distROI
%     The minimum distance in 8-connexity terms between the 'foreground'
%     element and the edges of the volume.
%   d8_matrix
%     The matrix of the 8-distances in the matrix.

dimMatrix = size(matrix);

%% Initialisation of the matrix of distances
% the background is assumed to have the value in the corner first voxel
value_bkgnd = matrix(1, 1, 1);
d8_matrix = zeros(dimMatrix);
d8_matrix(matrix == value_bkgnd) = inf;
d8_matrix = padarray(d8_matrix, [1, 1, 1], inf, 'post');
d8_matrix = padarray(d8_matrix, [1, 1, 1], inf, 'pre');

%% Scans
tic
% Span from upper left 
for i = 2:dimMatrix(1)+2
    for j = 2:dimMatrix(2)+1
        for k = 2:dimMatrix(3)+1
            tmp = [d8_matrix(i, j, k)  min(d8_matrix(i-1, j-1:j+1, k-1:k+1) + 1, [], 'all')  min(d8_matrix(i, j-1, k-1:k+1) + 1, [], 'all')  d8_matrix(i, j, k-1) + 1];
            d8_matrix(i, j, k) = min(tmp);
        end
    end
end

% Span from lower right
for i = dimMatrix(1)+1:-1:1
    for j = dimMatrix(2)+1:-1:2
        for k = dimMatrix(3)+1:-1:2
            tmp = [d8_matrix(i, j, k)  min(d8_matrix(i+1, j-1:j+1, k-1:k+1) + 1, [], 'all')  min(d8_matrix(i, j+1, k-1:k+1) + 1, [], 'all')  d8_matrix(i, j, k+1) + 1];
            d8_matrix(i, j, k) = min(tmp);
        end
    end
end
toc

% Truncate the distance matrix to fit the size of the initial matrix
d8_matrix = d8_matrix(2:dimMatrix(1) + 1, 2:dimMatrix(2) + 1, 2:dimMatrix(3) + 1);

%% Determination of the lowest distance on the edges
[distROI, ind] = min([min(d8_matrix(1, :, :), [], 'all') ...
            min(d8_matrix(:, 1, :), [], 'all') ...
            min(d8_matrix(:, :, 1), [], 'all') ...
            min(d8_matrix(:, :, end), [], 'all') ...
            min(d8_matrix(:, end, :), [], 'all') ...
            min(d8_matrix(end, :, :), [], 'all')]);
disp('end.')
end