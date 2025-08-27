%% pinv_recon_time_estimate_main
% uses script to estimate time it will take to do pinv-recon
% uses the function get_recon_time_estimate.m
% estimate is specific to device hardware
% this script also does the actual pinv-recon, to compare actual and
% theoretical results
%
% Sofia Pearson
% 01/07/2025

%% Example code to call:
% RUN FROM FOLDER Pinv-Recon_Dev;

clear all
clc

% Get folder Pinv-Recon_Dev and add all subfolders to path:
currentFolder = fileparts(mfilename('fullpath'));
targetFolder = 'Pinv-Recon_Dev';
while true
    [parentFolder, currentName] = fileparts(currentFolder);

    if strcmp(currentName, targetFolder)
        % Found
        rootDirectory = currentFolder;
        % Add Pinv-Recon_Dev and subfolders to the path
        addpath(genpath(rootDirectory));   
        % Go into Pinv-Recon_Dev folder
        cd(rootDirectory)
        break;
    elseif isempty(parentFolder) || strcmp(currentFolder, parentFolder)
        % Not found
        error('Pinv-Recon_dev Folder not found in parent hierarchy. Startup condition failed. Closing app.');        
    else
        % Move one level up
        currentFolder = parentFolder;
    end
end

%method
mode = 'cholesky';

%load data
wfn = fullfile(pwd, "Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/spiral_1h_fov240_mtx64_arms4_kdt4_gmax19_smax119_dur6p1_blncd.mat");
wfn = char(wfn);
data = load(fullfile(pwd, "Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/dd.mat"));
data = data.dd;
b0 = load(fullfile(pwd,"Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/fieldmap.mat"));
b0 = b0.fieldmap;
%sens
sens = load(fullfile(pwd, 'Pinv-Recon_App/Data_for_app/ScannerData(Real)/pinv_b0_sens_data/coil_sense_map.mat'));
sens = sens.rel_coil_sense;

%recon time, no correction maps
[recon_time_estimate_no_corrections] = get_recon_time_estimate(mode, 0, 'laptop', wfn, data, 0, 0 );
%recon time, with correction maps
%[recon_time_estimate_yes_corrections] = get_recon_time_estimate(mode, 0, 'laptop', wfn, data, 1, 1 );

%Do actual pinv!

tic,
[~,bbabs]=pinv_recon(data,wfn, 'mode', mode);
actual_recon_time_no_maps = toc;

% tic,
% [~,bbabs_with_maps]=pinv_recon(data,wfn,'b0', b0, 'sens', sens, 'mode', mode);
% actual_recon_time_with_maps = toc;

%results:
sprintf('No maps: Estimate time = %d seconds. Actual time = %d seconds.', recon_time_estimate_no_corrections, actual_recon_time_no_maps),
%sprintf('Yes maps: Estimate time = %d seconds. Actual time = %d seconds.', recon_time_estimate_yes_corrections, actual_recon_time_with_maps),

bbabs = bbabs(:,:,:,1);
figure, mat2montage(bbabs)