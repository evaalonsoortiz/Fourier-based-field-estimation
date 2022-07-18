function sus_lowRes = sub_sample(sus)
n_plots = 4;
sigma = 0.5;
sus_section = squeeze(sus(:, :, 65));
figure;
subplot(n_plots, 2, 1); imagesc(sus_section); colorbar; title('true susceptibility map');
%% Low-pass filtering
averaging = 1/9 * [1 1 1; 1 1 1; 1 1 1];
sus_section_averaging = conv2(sus_section, averaging, 'same');
subplot(n_plots, 2, 3); imagesc(sus_section_averaging); colorbar; title('averaging filtered susceptibility map');

sus_section_gaussian = imgaussfilt(sus_section, sigma, 'padding', 'symmetric');
subplot(n_plots, 2, 5); imagesc(sus_section_gaussian); colorbar; title('gaussian filtered susceptibility map');

sus_section_median = medfilt2(sus_section, [3, 3]);
subplot(n_plots, 2, 7); imagesc(sus_section_median); colorbar; title('gaussian filtered susceptibility map');

%% Down sampling

sus_lowRes = sus_section(1:2:end, 1:2:end);
subplot(n_plots, 2, 2); imagesc(sus_lowRes); colorbar; title('down-sampled susceptibility map');
fprintf('normal : %u \n',length(unique(sus_lowRes)))
disp(unique(sus_lowRes))


sus_lowRes = sus_section_averaging(1:2:end, 1:2:end);
subplot(n_plots, 2, 4); imagesc(sus_lowRes); colorbar; title('down-sampled averaging filtered susceptibility map');
fprintf('averaging : %u \n',length(unique(sus_lowRes)))

sus_lowRes = sus_section_gaussian(1:2:end, 1:2:end);
subplot(n_plots, 2, 6); imagesc(sus_lowRes); colorbar; title(sprintf('down-sampled gaussian filtered susceptibility map, sigma : %0.2f', sigma));
fprintf('gaussian sigma %0.2f : %u \n', sigma, length(unique(sus_lowRes)))

sus_lowRes = sus_section_median(1:2:end, 1:2:end);
subplot(n_plots, 2, 8); imagesc(sus_lowRes); colorbar; title('down-sampled median filtered susceptibility map');
fprintf('median : %u \n',length(unique(sus_lowRes)))
disp(unique(sus_lowRes))



figure;
sus_lowRes = sus_section(1:2:end, 1:2:end);
subplot(n_plots, 2, 1); imagesc(sus_lowRes); colorbar; title('down-sampled susceptibility map');

sus_lowRes = sus_section_averaging(1:2:end, 1:2:end);
subplot(n_plots, 2, 3); imagesc(sus_lowRes); colorbar; title('down-sampled averaging filtered susceptibility map');

sus_lowRes = sus_section_gaussian(1:2:end, 1:2:end);
subplot(n_plots, 2, 5); imagesc(sus_lowRes); colorbar; title('down-sampled gaussian filtered susceptibility map');

sus_lowRes = sus_section_median(1:2:end, 1:2:end);
subplot(n_plots, 2, 7); imagesc(sus_lowRes); colorbar; title('down-sampled median filtered susceptibility map');

sus_lowRes = sus_section(1:2:end, 1:2:end);
subplot(n_plots, 2, 2); imagesc(sus_lowRes(57:75, 102:111)); colorbar; title('down-sampled susceptibility map');

sus_lowRes = sus_section_averaging(1:2:end, 1:2:end);
subplot(n_plots, 2, 4); imagesc(sus_lowRes(57:75, 102:111)); colorbar; title('down-sampled averaging filtered susceptibility map');

sus_lowRes = sus_section_gaussian(1:2:end, 1:2:end);
subplot(n_plots, 2, 6); imagesc(sus_lowRes(57:75, 102:111)); colorbar; title('down-sampled gaussian (sigma 0.3) filtered susceptibility map');

sus_lowRes = sus_section_median(1:2:end, 1:2:end);
subplot(n_plots, 2, 8); imagesc(sus_lowRes(57:75, 102:111)); colorbar; title('down-sampled median filtered susceptibility map');



end