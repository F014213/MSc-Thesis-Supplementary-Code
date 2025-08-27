%% b0_structural_phantom_RESULTS_ANALYSIS

% Sofia Pearson
% 04/08/2025

%==========================================================================
%% LOAD DATA
% Load wavefrom, b0 map, and results
wfn = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/b0_example/waveform.mat');
wf = load_waveform(wfn);
b0map_name= fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/b0_example/fieldmap.mat');

load(b0map_name)
load(fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/b0_example/results.mat'));
%==========================================================================

% Plot k-space trajectory:
ktraj_plot = figure;
set(ktraj_plot, 'Units', 'inches', 'Position', [1, 1, 6.9, 4]); % [left bottom width height]

ktraj_plot = tiledlayout(1,2);
nexttile
plot(wf.k(:,:,1), wf.k(:,:,2)), axis square,
ax = gca;
set(ax, 'TickDir', 'out', 'TickLength', [0.02 0.02]); clear ax;
xlabel('kx', 'FontSize', 10), ylabel('ky', 'FontSize', 10), title('k-space Trajectory', 'FontSize', 12);

nexttile
imagesc(fieldmap), colormap gray, title('Off-Resonance Map', 'FontSize', 12);
axis square, axis off;

% save plot
% filename = 'C:\Users\sofia\OneDrive\Desktop\ktraj_b0map_dataset1.pdf';
% exportgraphics(ktraj_plot, filename,'Resolution', 300);



%% -------------------------------
% x = 0:1:(length(wf.t)-1);
% figure, plot(x, wf.t);

spiral_time = wf.t_arm;
num_arms = wf.nviews;
%=================================
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
%% coil combine
if ncoils>1
    IMG_wMaps = sqrt(mean(IMG_wMaps.*conj(IMG_wMaps),5));
end


%==========================================================================

%% Two Plots
recon_plots = figure;
colormap gray;
set(recon_plots, 'Units', 'inches', 'Position', [1, 1, 6.9, 4]); % [left bottom width height]

recon_plots = tiledlayout(1,2,'TileSpacing', 'compact', 'Padding', 'compact');
nexttile
imagesc(IMG), title('Recon', 'FontSize', 12), axis square, axis off;

nexttile
imagesc(IMG_wMaps), title('Recon + Off-Res', 'FontSize', 12), axis square, axis off;
hold off;


% save plot
% filename = 'C:\Users\sofia\OneDrive\Desktop\all_recons_dataset4.pdf';
% exportgraphics(recon_plots, filename,'Resolution', 300);
% 

%% ========================================================================
%% IMAGE SHARPNESSSSSS:

%imgradient3

[Gmag, ~, ~] = imgradient3(IMG);
[Gmag_wMaps, ~, ~] = imgradient3(IMG_wMaps);

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


%% Fuse Images? - DO THIS LATER
%fuse IMG and fieldmap:

%background image
bg = IMG;
bg = mat2gray(bg);  % scale to [0,1]

B0_map = imresize(fieldmap,[96 96]);
B0_map = B0_map/(max(abs(fieldmap), [], 'all'));

% Define center value
centerVal = 0;
range = max(abs(B0_map(:)-centerVal));

figure;
imshow(bg, []);  
hold on;

h = imagesc(B0_map, [centerVal - range, centerVal + range]);

% Build custom blue-white-red colormap
n = 256;  % number of colors
half = n/2;

% From blue (0,0,1) to white (1,1,1)
blue_to_white = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];

% From white (1,1,1) to red (1,0,0)
white_to_red = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];

% Concatenate to form full colormap
cmap = [blue_to_white; white_to_red];

% Apply colormap
colormap(cmap);

% %% Now fuse the images! 
 set(h, 'AlphaData', 0.7 * (B0_map ~= centerVal));  


% alphaMask = abs(B0_map - centerVal) / range; 
% set(h, 'AlphaData', alphaMask);

%% =========

%% SRF
SRF = PinvReconData.SRFMatrixGradientOnly;
SRF_b0 = PinvReconData.SRFMatrix;

SRF_mean = abs(mean(SRF, 'all'));
SRF_b0_mean = abs(mean(SRF_b0, 'all'));

SRF_glob_max = max([max(abs(SRF)), max(abs(SRF_b0))]);
SRF_glob_min = min([min(abs(SRF)), min(abs(SRF_b0))]);

% Plot k-space trajectory:
SRF_maps = figure;
set(SRF_maps, 'Units', 'inches', 'Position', [1, 1, 6.9, 4]); % [left bottom width height]

SRF_maps = tiledlayout(1,2);
nexttile
imagesc(abs(SRF)), axis square, axis off; colorbar;
title('SRF Map Gradient-Only', 'FontSize', 12);
clim([SRF_glob_min, SRF_glob_max]);

nexttile
imagesc(abs(SRF_b0)), axis square, axis off; colormap hot; colorbar;
title('SRF Map with Off-Res', 'FontSize', 12);
clim([SRF_glob_min, SRF_glob_max]);

% save plot
% filename = 'C:\Users\sofia\OneDrive\Desktop\ktraj_b0map_dataset1.pdf';
% exportgraphics(ktraj_plot, filename,'Resolution', 300);
