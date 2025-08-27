%% pinv_xenon_b0_lungs_RESULTS_ANALYSIS

% Sofia Pearson
% 04/08/2025

%==========================================================================
%% LOAD DATA
% Load wavefrom, b0 map, and results
wfn = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_xenon_b0_data/15mSpiral1ArmXenon.mat');
wf = load_waveform(wfn);
b0map_name= fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_xenon_b0_data/B0map.mat');

load(b0map_name)
load(fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_xenon_b0_data/results.mat'));
%==========================================================================

b0map_flat = mat2montage_for_app(B0);

% Plot B0 + ktraj 

map_plots = figure;
colormap gray;
set(map_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

map_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(wf.k(:,:,1), wf.k(:,:,2)), axis square,
ax = gca;
set(ax, 'TickDir', 'out', 'TickLength', [0.02 0.02]); clear ax;
xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('k-space Trajectory', 'FontSize', 12);

nexttile
imagesc(b0map_flat), title('Off-Resonance Map', 'FontSize', 12), axis equal, axis tight, axis off;


filename = 'C:\Users\sofia\OneDrive\Desktop\ktraj_b0_xenon.pdf';
exportgraphics(map_plots, filename,'Resolution', 300);


%% -------------------------------

% x = 0:1:(length(wf.t)-1);
% figure, plot(x, wf.t);

spiral_time = wf.t_arm;
num_arms = wf.nviews;


%% -------------------------------

%==========================================================================
% analyse results:

IMG = PinvReconData.ImageGradientOnly;
IMG_wMaps = PinvReconData.ImageWithMaps;

sz =size(IMG);    % get the dimensions of image (x,y,z,t,coils)
numDims = length(sz);
if numDims == 3                
    nslices = sz(3);
elseif numDims == 4             
    nslices = sz(3);        
    ntimesteps = sz(4);
elseif numDims == 5  
    nslices = sz(3);        
    ntimesteps = sz(4);
    ncoils = sz(5);
end          
if ~isempty(ntimesteps)
    disp('Averaging over the time steps.')
    IMG = sum(IMG,4); % sum the timesteps
end
%% coil combine
if ncoils>1
    IMG = sqrt(mean(IMG.*conj(IMG),5));
end


sz =size(IMG_wMaps);    % get the dimensions of image (x,y,z,t,coils)
numDims = length(sz);
if numDims == 3                
    nslices = sz(3);
elseif numDims == 4             
    nslices = sz(3);        
    ntimesteps = sz(4);
elseif numDims == 5  
    nslices = sz(3);        
    ntimesteps = sz(4);
    ncoils = sz(5);
end          
if ~isempty(ntimesteps)
    disp('Averaging over the time steps.')
    IMG_wMaps = sum(IMG_wMaps,4); % sum the timesteps
end


flat_IMG = mat2montage_for_app(IMG);
flat_IMG_wMaps = mat2montage_for_app(IMG_wMaps);


%% Plot of Recons
recon_plots = figure;
colormap gray;
set(recon_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

recon_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(flat_IMG), title('Recon', 'FontSize', 12), axis equal, axis tight, axis off;

nexttile
imagesc(flat_IMG_wMaps), title('Recon + Off-Res', 'FontSize', 12), axis equal, axis tight, axis off;


filename = 'C:\Users\sofia\OneDrive\Desktop\all_recons_xenon.pdf';
exportgraphics(recon_plots, filename,'Resolution', 300);

%%=========================================================================4%


% Analysis of Slice 7:
IMG_7 = abs(IMG(:,:,7));
IMG_wMaps_7 = abs(IMG_wMaps(:,:,7));

global_min_7 = min([ min(IMG_7(:)), min(IMG_wMaps_7(:)) ]);
global_max_7 = max([max(IMG_7(:)), max(IMG_wMaps_7(:)) ]);

slice_7_plots = figure;
colormap gray;
set(slice_7_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]

slice_7_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(IMG_7), title('Slice 7', 'FontSize', 12), axis equal, axis tight, axis off; clim([global_min_7, global_max_7]);

nexttile
imagesc(IMG_wMaps_7), title('Slice 7 + Off-Res', 'FontSize', 12), axis equal, axis tight, axis off; clim([global_min_7, global_max_7]);


filename = 'C:\Users\sofia\OneDrive\Desktop\slice_7_xenon.pdf';
exportgraphics(slice_7_plots, filename,'Resolution', 300);

%b0 map for slice 7:
B0_s7 = B0(:,:,7);

B0_s7_plot = figure;
colormap gray;
set(B0_s7_plot, 'Units', 'inches', 'Position', [1, 1, 3.5, 3.5]); % [left bottom width height]
imagesc(B0_s7), title('Slice 7 Off-Resonance Map', 'FontSize', 12), axis equal, axis tight, axis off; 


filename = 'C:\Users\sofia\OneDrive\Desktop\B0_s7_xenon.pdf';
exportgraphics(B0_s7_plot, filename,'Resolution', 300);


%% imgradiet3 of image

[Gmag, ~, ~] = imgradient3(abs(IMG));
[Gmag_wMaps, ~, ~] = imgradient3(abs(IMG_wMaps));

mean_grad = mean(Gmag, "all");
mean_grad_wMaps = mean(Gmag_wMaps, "all");

%plot gradients
glob_max_2d = max([max(Gmag), max(Gmag_wMaps)]);
glob_min_2d = min([min(Gmag), min(Gmag_wMaps)]);

figure,
tiledlayout(1,2)
nexttile
imagesc(Gmag), axis square, axis off; colorbar;
title('Gradient magnitude')
clim([glob_min_2d, glob_max_2d]);

nexttile
imagesc(Gmag_wMaps), axis square, axis off; colorbar;
title('Gradient magnitude')
clim([glob_min_2d, glob_max_2d]);





%% ========================================================================
% % Analysis of Slice 6:
% IMG_6 = abs(IMG(:,:,6));
% IMG_wMaps_6 = abs(IMG_wMaps(:,:,6));
% 
% global_min_6 = min([ min(IMG_6(:)), min(IMG_wMaps_6(:)) ]);
% global_max_6 = max([max(IMG_6(:)), max(IMG_wMaps_6(:)) ]);
% 
% slice_6_plots = figure;
% colormap gray;
% set(slice_6_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]
% 
% slice_6_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
% nexttile
% imagesc(IMG_6), title('Slice 6', 'FontSize', 12), axis equal, axis off; clim([global_min_6, global_max_6]);
% 
% nexttile
% imagesc(IMG_wMaps_6), title('Slice 6 + Off-Res', 'FontSize', 12), axis equal, axis off; clim([global_min_6, global_max_6]);
% 
% sub = IMG_6 - IMG_wMaps_6;
% figure, imagesc(sub), axis equal, colormap gray; colorbar;
% 
% 
% %slice 3 coil sense and slice 3 B0 map,
% slice_3_b0_map = fieldmap(:,:,3);
% coil_sense_slice_3 = rel_coil_sense(:,:,3,:);
% coil_sense_slice_3 = squeeze(coil_sense_slice_3);
% flat_coil_sense_slice_3 = mat2montage_for_app(coil_sense_slice_3);
% 
% slice_3_b0_plot = figure;
% 
% colormap gray;
% set(slice_3_b0_plot, 'Units', 'inches', 'Position', [1, 1, 6.9, 5]); % [left bottom width height]
% 
% slice_7_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
% nexttile
% imagesc(slice_3_b0_map), title('Slice 3 Off-Resonance Map', 'FontSize', 12), axis square, axis off;
% 
% nexttile
% imagesc(flat_coil_sense_slice_3 ), title('Slice 3 Coil Sensitivity Maps', 'FontSize', 12), axis square, axis off;
% hold off;
% 
% 
% filename = 'C:\Users\sofia\OneDrive\Desktop\slice_3_maps_dataset3.pdf';
% exportgraphics(slice_7_plots, filename,'Resolution', 300);



%% SRF
SRF = PinvReconData.SRFMatrixGradientOnly;
SRF_b0_coil = PinvReconData.SRFMatrix;

SRF_b0_coil_slice7 = SRF_b0_coil(:,:,3);

SRF_mean = abs(mean(SRF, 'all'));
SRF_b0_coil_mean = abs(mean(SRF_b0_coil, 'all'));

% SRF_glob_max = max([max(abs(SRF)), max(abs(SRF_b0_coil_slice7))]);
% SRF_glob_min = min([min(abs(SRF)), min(abs(SRF_b0_coil_slice7))]);

% Plot k-space trajectory:
SRF_maps = figure;
set(SRF_maps, 'Units', 'inches', 'Position', [1, 1, 6.9, 4]); % [left bottom width height]

SRF_maps = tiledlayout(1,2);
nexttile
imagesc(abs(SRF)), axis square, axis off; colorbar;
title('SRF Map Gradient-Only', 'FontSize', 12);
% clim([SRF_glob_min, SRF_glob_max]);

nexttile
imagesc(abs(SRF_b0_coil_slice7)), axis square, axis off; colormap hot; colorbar;
title('SRF Map with Off-Res', 'FontSize', 12);
% clim([SRF_glob_min, SRF_glob_max]);

% save plot
% filename = 'C:\Users\sofia\OneDrive\Desktop\ktraj_b0map_dataset1.pdf';
% exportgraphics(ktraj_plot, filename,'Resolution', 300);
