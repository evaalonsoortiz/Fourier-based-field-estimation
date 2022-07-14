%% Script to test subsample and the effect on resolution on the Zubal phantom


% generate susceptibility distribution for the modified Zubal phantom
zubal_sus_dist = Zubal('zubal_EAO.nii');

% Properties of the phantom : dimensions 256x256x128
% Resolution 1.1x1.1x1.4