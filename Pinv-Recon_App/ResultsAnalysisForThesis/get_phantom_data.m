function [phantom_data, b0Map] = get_phantom_data(wfn, varargin)
%        phantom_data = get_phantom_data(wfn, 'npix', npix, 'b0', 1 or 0, 'save', 1 or 0)

%% Prepare synthetic MR data for Pinv-Recon using custom waveform and Shepp-Logan Head Phantom
% 
% Inputs       
%             wfn   waveform name     
%
% Optional inputs
%            npix   reconstructed image matrix size                    ([])
%              b0   also get synthetic b0 map and apply b0              (0)
%                   inhomogeneity effects to the data. 0 no b0, 1 yes b0
%            save   set to 1 or 0. Saves the synthetic MR data to the   (0) 
%                   same folder as the waveform, as a .mat file.
% Outputs
%    phantom_data   synthetic MR data (k-space data) based on Shepp-Logan
%                   head pahntom and the input waveform. 
%
% 07/2025   Sofia Pearson

%% Example function call: (run from folder Pinv-Recon_Dev)
% wfn = fullfile(pwd, 'Pinv-Recon_App/Data_for_app/TrajectoryTestingData(Phantom)/n_spiral_arms/01arm/spiral_1arm.mat');
% [phantom_data, b0Map] = get_phantom_data(wfn, 'npix', 128, 'b0', 1, 'save', 1);


    % GET INPUTS:
    % ---------------------------------------------------------------------
    
    % Set default values
    npix = [];
    saveFlag = 0; %don't save
    b0 = 0; % don't apply b0 effects
    b0Map = [];   

    % Loop through varargin
    i = 1;
    while i <= length(varargin)
        if ischar(varargin{i}) || isstring(varargin{i})
            switch lower(varargin{i})
                case 'npix'
                    npix = varargin{i+1};
                    i = i + 2;
                case 'b0'
                    b0 = varargin{i+1};
                    i = i + 2;
                case 'save'
                    saveFlag = varargin{i+1};
                    i = i + 2;
                otherwise
                    error('Unknown parameter: %s', varargin{i});
            end
        else
            error('Parameter names must be strings.');
        end
    end

    % ---------------------------------------------------------------------
    % GET DATA:

    % Step 1: Load Waveform:
    wf = load_waveform(wfn);
    k = wf.k;
    
    % if npix is not provided:
    if isempty(npix)
        if ~isfield(wf,'mtx') && ~isfield(wf,'npix')
            error('Acquisition matrix size is required in wf.mtx or wf.npix');
        end
        try mtx_acq=wf.mtx; catch; mtx_acq=wf.npix; end
        if length(mtx_acq)==1
            mtx_acq=[mtx_acq mtx_acq];
        end
        npix=mtx_acq;
    end

    % Step 2: Generate Phantom:
    P = phantom(npix);
    P(P>0.75)=0.75; P=single(P/max(P(:))); %normalise

    %Step 3: Get Encode Matrix:
    
%%  %image coords (Expand to 1D and 3D coords in future)
    [XM,YM]=ndgrid(single(-(npix/2-1):npix/2));
    
    % scale from +-0.5 to +-pi
    Kx=single(k(:,:,1))./0.5*pi./npix.*npix;
    Ky=single(k(:,:,2))./0.5*pi./npix.*npix;
    
     
    % Encoding matrix exp(i*k*r)
    E=exp(1i*(Kx(:)*(XM(:)).'+Ky(:)*(YM(:)).'));

    %Step 4: ADD B0 MAP (optiona;)
    if b0 == 1
        %get b0 map
        time = wf.t;
        xi = linspace(-1,1,npix+1);
        xi = xi(1:end-1);
        [~,Y_b0]=ndgrid(xi, xi);
        b0Map = 125*Y_b0.^2-30;
        %encode matrix for b0 map
        E_B0 = exp(-1i*2*pi*time.'*(b0Map(:).')); 
        %update overall encode matrix
        E = E.*E_B0;
    end

    % Step 5: Get Data:
    phantom_data = E * P(:);
    phantom_data = phantom_data.';

    % Step 6: Save Data:

    if saveFlag == 1
        [wfn_folder, wfn_name] = fileparts(wfn);   

        if b0 == 0
            fname = [wfn_folder, '/data_for_', wfn_name, '.mat'];
        elseif b0 == 1
            fname = [wfn_folder, '/data_with_b0_for_', wfn_name, '.mat'];
        end

        fname= char(fname);       
    
        if exist(fname,'file')
            fprintf('Phantom data file for this waveform already exists as: \n %s', fname)
        else
            save(fname, "phantom_data")
            fprintf('Phantom data file for this waveform successfully saved as: \n %s', fname)
        end
        % Save b0 Map - if present
        if b0 == 1
            b0_fname = [wfn_folder, '/b0_map_for_', wfn_name, '.mat'];
            b0_fname= char(b0_fname);
            if exist(b0_fname,'file')
                fprintf('\n B0 Map file for this waveform already exists as: \n %s', fname)
            else
                save(b0_fname, "b0Map")
                fprintf('\nB0 Map file for this waveform successfully saved as: \n %s', fname)
            end
        end
    end
end