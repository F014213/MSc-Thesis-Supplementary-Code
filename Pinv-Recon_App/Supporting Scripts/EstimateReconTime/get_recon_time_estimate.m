%% Pinv-Recon time estimator:
% 
% Inputs - pinv method - 'cholesky', 'svd', 'eig', 'qr'
%        - useGPU - 1 or 0
%        - mtx_reco - size of reconstruction matrix
%        - waveform size ( length(kx))
%        - nslices
%        - hardware - 'laptop' or 'workstation'
%        - b0 is just a yes or no, 1 or 0
%        - coilsense is just a yes or no, 1 or 0
%
% Output - Recon Time Estimate - with and without correction
%
% 27/06/2025
% Sofia Pearson

%% To run the function:
% see script pinv_recon_time_estimate_main

function [recon_time_estimate] = get_recon_time_estimate(method, useGPU, hardware, wfn, data, coilsense, b0)


%% code from pinv_recon:---------------------------------------------------

%% reading waveform file
if isstruct(wfn)
    wf = wfn;
    wfn='unknown';Fdir='unknown';name='unknown';
    warning('Waveform provided in struct format. Will not be attempting to read or write reconstruction matrix.')
else
        % load waveform structure from fidall.e
    [Fdir,name,ext]=fileparts(wfn);
    wfn=[Fdir '/' name '.mat'];
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    wf = load_waveform(wfn,true,true);       % load waveform
    readfile=1;
end

%% ensure mtx_reco is a 1x3
[ncoils,npts,ntimesteps,nslices,nspec]=size(data);
try mtx_acq=wf.mtx; catch; mtx_acq=wf.npix; end
if length(mtx_acq)==1;mtx_acq=[mtx_acq mtx_acq];end
mtx_reco=mtx_acq;
if length(mtx_acq)==1; mtx_acq=[mtx_acq mtx_acq nslices];end
if length(mtx_acq)==2; mtx_acq=[mtx_acq nslices];end
if length(mtx_reco)==1; mtx_reco=[mtx_reco mtx_reco nslices];end
if length(mtx_reco)==2; mtx_reco=[mtx_reco nslices];end

%% get reconstruction matrix
k=wf.k;
if ~isreal(k)
    if size(k,1)~=1, error('size(k,1)~=1'); end
    k_temp=k;
    k(1,:,1) = real(k_temp);
    k(1,:,2) = imag(k_temp);
    clear k_temp
end

%%-------------------------------------------------------------------------


% calibration 
    switch hardware
        case 'laptop'
            Ni = [16, 32, 48, 64];    
        case 'workstation'
            Ni = [48, 64, 96];
    end
    
    if useGPU == 1
        gpu = 1;
    else 
        gpu = 0; % assume no gpu
    end
    

    switch method

        %% ADD IN GPU CASES - AND FIX GPU ISSUE
        case 'svd'
            for i1=1:length(Ni)
                N=Ni(i1)^2; 
                E=randn(N,N,'single')+1i*randn(N,N,'single'); 
                % CPU
                if Ni(i1)<257
                    tic, [U,S,V]=svd(E,'econ'); inversion_times(i1)=toc, clear U S V;                     
                end                           
            end

        case 'qr'
            for i1=1:length(Ni)
                N=Ni(i1)^2; 
                E=randn(N,N,'single')+1i*randn(N,N,'single'); 
                % CPU
                if Ni(i1)<257
                    tic, [Q,R,p1]=qr(E,'econ','vector'); inversion_times(i1)=toc, clear Q R p1;                    
                end                           
            end

        case 'eig'
            for i1=1:length(Ni)
                N=Ni(i1)^2; 
                E=randn(N,N,'single')+1i*randn(N,N,'single'); 
                % CPU
                tic, EHE=E'*E+1.0*eye(N,N); 
                [V,D]=eig(EHE); inversion_times(i1)=toc, clear V D; 
                
                                                 
            end
        case 'cholesky'
            for i1=1:length(Ni)
                N=Ni(i1)^2; 
                E=randn(N,N,'single')+1i*randn(N,N,'single'); 
                % CPU
                tic, EHE=E'*E+1.0*eye(N,N);         
                L=chol(EHE); inversion_times(i1)=toc, clear L;                           
            end

    end

    % actual size of E?
    kx = k(1,:,1);
    nRows = length(kx);

    if coilsense ==1
        nRows = nRows * ncoils; %only multiply by ncoils if we actually use the coil sensititvity maps? 
    end

    nCols = mtx_reco(1).^2; % this wont work, fix it
    
    sizes = Ni; %.^2;
    coeffs = polyfit(sizes, inversion_times, 3);  % Fit a cubic polynomial
    % % find number of elements in E:
    % n_elements = nRows*nCols;
    % recon_time_estimate = polyval(coeffs, n_elements);
    
    n_elements = round(sqrt(nRows*nCols));
    recon_time_estimate = polyval(coeffs, n_elements);

    %repeat if b0 and coil sense, has to be done per slice
    if b0 == 1 || coilsense ==1
        recon_time_estimate = recon_time_estimate * nslices;
    end
            
    % for i1=1:length(Ni)
    %     N=Ni(i1)^2; 
    %     E=randn(N,N,'single')+1i*randn(N,N,'single'); 
    %     % CPU
    %     if Ni(i1)<257
    %         tic, [U,S,V]=svd(E,'econ'); cpuSVD(i1)=toc, clear U S V; 
    %         tic, [Q,R,p1]=qr(E,'econ','vector'); cpuQRD(i1)=toc, clear Q R p1;
    %     end
    % 
    %     tic, EHE=E'*E+1.0*eye(N,N); cpuEHE(i1)=toc, 
    %     tic, [V,D]=eig(EHE); cpuEIG(i1)=toc, clear V D; 
    %     tic, L=chol(EHE); cpuCHOL(i1)=toc, clear L;
    %     % GPU
    %     if Ni(i1)<193
    %         E=gpuArray(E);
    %         tic, [U,S,V]=svd(E,'econ'); gpuSVD(i1)=toc, clear U S V; 
    %         tic, [Q,R,p1]=qr(E,'econ','vector'); gpuQRD(i1)=toc, clear Q R p1;
    % 
    %         tic, EHE=E'*E+1.0*eye(N,N); gpuEHE(i1)=toc, clear E;
    %         tic, [V,D]=eig(EHE); gpuEIG(i1)=toc, clear V D; 
    %         tic, L=chol(EHE); gpuCHOL(i1)=toc, clear L;
    %     end
    % end
end


