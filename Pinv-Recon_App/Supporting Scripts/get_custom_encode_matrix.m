% Get encode matrix for Pinv-Recon_App
% 23/05/2025
% Sofia Pearson

% Making a function. Adapt for app later. 

% Input Vars---------------------------------------------------------------

%ntx=64;%64,96,128

%TO RUN:
% [flags, recon_img] = get_custom_encode_matrix(64, 'gradNonLin');

%--------------------------------------------------------------------------
function [flags, recon_img] = get_custom_encode_matrix(varargin)

    % addpath:
    addpath(genpath(pwd))

    CondNumb=30; 
    ntx = 64; %default value
    theta = 15; %degrees

    %find ntx in varargin (if present):
    for i = 1:length(varargin)
        if isnumeric(varargin{i})
            ntx = varargin{i}; 
            break;
        end
    end
    
    %----------------------------------------------------------------------
    % Generate k-space Trajectory:
    fov  =180; arms = 1; ksamp = 4; gmax = 33; smax = 120; nucleus = '1H'; acq_round2n = true; balanced = true;
    %from -0.5 to 0.5
    [i,~,t,~,~]=design_spiral_for_app(fov,ntx,arms,ksamp,gmax,smax,nucleus,acq_round2n,balanced);
    time = t;
    Kx = real(i);
    Ky = imag(i);
    %Check limits of k-space - ensure this is -1 to 1
    kmax = max([max(Kx), abs(min(Kx)), max(Ky), abs(min(Ky))]);
    %rescale trajectory to [-1,1]:
    Kx = Kx./kmax;
    Ky = Ky./kmax; 
    % kspace = [Kx' Ky'];
    
    %----------------------------------------------------------------------
    % Shepp-Logan phantom and image coords
    P = phantom('Modified Shepp-Logan',ntx); 
    x = linspace(-ntx,ntx,ntx);
    [X, Y] = ndgrid(x, x);
   
    %get gradient encode
    grad_encode = single(exp(1i*2*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).')));
    
    %----------------------------------------------------------------------
    % Get selected encode effects:

    %determine what we are encoding for:
    validStrings = {'b0', 'gradNonLin', 'coilSense','chemShift', 'rotation', 'offCtr'};
  
    % initialize flags
    flags = struct();
    for i = 1:length(validStrings)
        flags.(validStrings{i}) = false;
    end

    % Loop through inputs and check for matches (case-insensitive)
    for i = 1:length(varargin)
        if ischar(varargin{i}) || isstring(varargin{i})
            str = char(varargin{i});
            validatedStr = validatestring(str, validStrings, mfilename, 'input', i);
            flags.(validatedStr) = true;
        end
    end
    

    % Store selected encode matrices in a cell array
    encode_matrices = {grad_encode};

    %----------------------------------------------------------------------
   
    % Get encode matrices
        % TO INCLUDE OFFCTR AND ROTATION, MODIFY X AND Y BEFORE CALCULATING
    % OTHER EFFECTS??

     if flags.rotation
        rot_encode = imrotate(encode_matrices{1}, theta, 'bilinear', 'crop'); 
        %replace grad_encode with rot_encode
        encode_matrices{1} = rot_encode;
               
    end

    if flags.offCtr
        offCtr_encode = []; % ADD CODE
        encode_matrices{1} = offCtr_encode;
    end

    %-------------------------------------------------------------
    if flags.b0
        % get b0 map:
        Ynorm = Y./max([max(Y) abs(min(Y))]) ;
        B0Map = 125*Ynorm.^2-30;
        b0_encode = single(exp(1i*time(:)*(B0Map(:)).')); % if units is in Hertz;
        % fprintf('b0 encode size: %d .', size(b0_encode,2));
        encode_matrices{end+1} = b0_encode;
    end
    
    if flags.gradNonLin
        % get gradient non-lineearity map
        Xnorm = X./max([max(X) abs(min(X))]) ;        
        gradNonLin_Map = -(0.15*(Xnorm.^3));
        figure, imagesc(gradNonLin_Map), title('map'), axis square, colormap gray;
        gradNonLin_encode = single(exp(2i*pi*(Kx(:)*(X(:)+1*gradNonLin_Map(:)).'+Ky(:)*(Y(:)+0*gradNonLin_Map(:)).')));
        % REPLACE X AND Y WITH SHIFTED OR ROTATED COORDS IF NECESSARY
        %max_val = max(imag(gradNonLin_encode));
        encode_matrices{1} =  gradNonLin_encode; %replaces grad_encode since it includes this information
      
    end

    if flags.coilSense
        coilSense_encode =[]; % ADD CODE
        encode_matrices{end+1} = coilSense_encode;
    end

    if flags.chemShift
        chemShift_encode = []; % ADD CODE
        encode_matrices{end+1} = chemShift_encode;
    end

  

    encode_matrix = encode_matrices{1};

    for i = 2:length(encode_matrices)
        encode_matrix  =   encode_matrix .* encode_matrices{i}; % .* or * ?
    end
    
    kspace_data = encode_matrix*P(:);
    %----------------------------------------------------------------------
    % Reconstruct!
    fprintf('>> Computing gradient only reconstruction ... ');
    tic, 
    [U1,S1,V1]=svd((grad_encode),'econ');
    imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
    grad_recon_mat =V1*invS1*U1';
    toc, 
    grad_recon_img = reshape(grad_recon_mat*kspace_data, [ntx, ntx]);
    
    
    fprintf('>> Computing complete reconstruction ... ');
    tic,
    [U1,S1,V1]=svd((encode_matrix),'econ'); 
    imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
    recon_mat =V1*invS1*U1';
    toc, 
    recon_img = reshape(recon_mat*(kspace_data), [ntx, ntx]);


    fprintf('>> Computing TEST reconstruction (Grad)... ');
    tic,
    gradNonLin_encode = encode_matrices{end};
    [U1,S1,V1]=svd((gradNonLin_encode),'econ'); 
    imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
    recon_mat_test =V1*invS1*U1';
    toc, 
    recon_img_test = reshape(grad_recon_mat*(gradNonLin_encode*P(:)), [ntx, ntx]); %the incomplete recon

    recon_img_test_1 = reshape(recon_mat_test*(gradNonLin_encode*P(:)), [ntx, ntx]); %the recon
        

    % plotting
    figure
    tiledlayout(1,4)
    
    nexttile
    imagesc(abs(grad_recon_img)),  title('Grad Only Recon'), axis square, colormap gray;
    
    nexttile
    imagesc(abs(recon_img)),  title(sprintf('Complete Recon')), axis square, colormap gray;
    
    nexttile    
    imagesc(abs(recon_img_test)),  title(sprintf('grad only recon test')), axis square, colormap gray;

    nexttile 
    imagesc(abs(recon_img_test_1)),  title(sprintf('grad only recon test')), axis square, colormap gray;

    
end

% %%% TEST SCRIPT:
% 
% 
% CondNumb=30; 
% ntx = 64; %default value
% theta = 15; %degrees
% 
% %find ntx in varargin (if present):
% for i = 1:length(varargin)
%     if isnumeric(varargin{i})
%         ntx = varargin{i}; 
%         break;
%     end
% end
% 
% %----------------------------------------------------------------------
% % Generate k-space Trajectory:
% addpath('C:\Users\sofia\OneDrive\Desktop\Oxford_Molecular_Imaging\Matlab\Fidall\')
% fov  =180; arms = 1; ksamp = 4; fname = false; gmax = 33; smax = 120; nucleus = '1H';
% acq_round2n = true; do_rot_file = false; balanced = true;
% %from -0.5 to 0.5
% [k,dcf,t,ind,out,grad]=design_spiral(fov,ntx,arms,ksamp,fname,gmax,smax,nucleus,acq_round2n,do_rot_file,balanced);
% 
% time = t;
% Kx = real(k);
% Ky = imag(k);
% %Check limits of k-space - ensure this is -1 to 1
% kmax = max([max(Kx), abs(min(Kx)), max(Ky), abs(min(Ky))]);
% %rescale trajectory:
% fprintf('HERE')
% Kx = Kx./kmax;
% Ky = Ky./kmax; 
% kspace = [Kx' Ky'];
% 
% %----------------------------------------------------------------------
% % Shepp-Logan phantom and image coords
% P = phantom('Modified Shepp-Logan',ntx); 
% x = linspace(-ntx,ntx,ntx);
% [X, Y] = ndgrid(x, x);
% 
% %get gradient encode
% grad_encode = single(exp(1i*2*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).')));
% 
% 
% 
% 
% % Store selected encode matrices in a cell array
% encode_matrices = {grad_encode};
% 
% 
% % get b0 map:
% B0Map = 125*Y.^2-30;
% b0_encode = single(exp(1i*time(:)*(B0Map(:)).')); % if units is in Hertz;
% 
% 
% kspace_data = grad_encode*P(:);
% encode_matrix  =   grad_encode*
% for i = 2:length(encode_matrices)
%     encode_matrix  =   encode_matrix .* encode_matrices{i}; % .* or * ?
% end
% 
% %----------------------------------------------------------------------
% % Reconstruct!
% fprintf('>> Computing gradient only reconstruction ... ');
% tic, 
% [U1,S1,V1]=svd((encode_matrices{1}),'econ');
% imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
% grad_recon_mat =V1*invS1*U1';
% toc, 
% grad_recon_img = reshape(grad_recon_mat*kspace_data, [ntx, ntx]);
% 
% 
% fprintf('>> Computing complete reconstruction ... ');
% tic,
% [U1,S1,V1]=svd((encode_matrix),'econ'); 
% imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
% recon_mat =V1*invS1*U1';
% toc, 
% recon_img = reshape(recon_mat*(kspace_data), [ntx, ntx]);
% 
% % plotting
% figure
% tiledlayout(1,2)
% 
% nexttile
% imagesc(abs(grad_recon_img)),  title('Grad Only Recon'), axis square, colormap gray;
% 
% nexttile
% imagesc(abs(recon_img)),  title(sprintf('Complete Recon')), axis square, colormap gray;
% 
% 
% 
% 
