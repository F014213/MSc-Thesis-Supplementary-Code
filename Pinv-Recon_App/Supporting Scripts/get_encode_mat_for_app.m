% Get encode matrix for Pinv-Recon_App
% 19/05/2025
% Sofia Pearson
% Inputs: Selected Encode matrices. 

%in app, initalise all encode matrices as empty arrays, that are updated
%when the boxes are ticked. 

%% TEST RUN: ==============================================================
clear variables variables
%==============================================
    fprintf('>> Simulating ...');
%==============================================

% define empty encode matrices
gradient_encode = [];
offCtr_encode = [];
b0_encode = [];
gradNonLin_encode = [];
coilSense_encode = [];
chemShift_encode = [];
rot_encode = [];
spatioTemporal_encode = [];

%other input vars:
CondNumb=30; 
ntx=64;%64,96,128

%% Generate Trajectory:
addpath('C:\Users\sofia\OneDrive\Desktop\Oxford_Molecular_Imaging\Matlab\Fidall\')
fov  =180;
npix = 128;
arms = 1;
ksamp = 4;
fname = false;
gmax = 33;
smax = 120;
nucleus = '1H';
acq_round2n = true;
do_rot_file = false;
balanced = true;
%from -0.5 to 0.5
[k,dcf,t,ind,out,grad]=design_spiral(fov,npix,arms,ksamp,fname,gmax,smax,nucleus,acq_round2n,do_rot_file,balanced);
Kx = real(k);
Ky = imag(k);
%Check limits of k-space - ensure this is -1 to 1
kmax = max([max(Kx), abs(min(Kx)), max(Ky), abs(min(Ky))]);
%rescale trajectory:
Kx = Kx./kmax;
Ky = Ky./kmax;


% Shepp-Logan phantom:
imgPhantom = phantom('Modified Shepp-Logan',ntx); 

% Generate the B0 maps: 
x = linspace(-ntx,ntx,ntx);
[X, Y] = ndgrid(x, x);
kspace = [Kx' Ky'];

B0Map=125*Y.^2-30; % why are these scale factors chosen? %(ASSUME Hz)
% field map range across 125 Hz total variation with a -30 Hz global shift at the center


% Gradient Encoding
gradient_encode = single(exp(1i*2*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).')));

% Off Resonance Encoding
%b0_encode = single(exp(1i*time(:)*(B0Map(:)).')); % if units is in Herts

%run function to get encode
Encode = get_custom_encode_matrix(gradient_encode, offCtr_encode, b0_encode, gradNonLin_encode, coilSense_encode, chemShift_encode, rot_encode, spatioTemporal_encode);
    
%get recon
tic, [U,S,V]=svd((Encode),'econ'); toc, 
diagS=diag(S); imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
tic, Recon =V*invS*U'; toc, 
IMG = reshape(Recon*(Encode*imgPhantom(:)), [ntx, ntx]);
figure, imagesc(abs(IMG)), title('Pinv Recon'), axis square, 
%% ========================================================================


%% FUNCTION
function ENCODE = get_custom_encode_matrix(gradient_encode, offCtr_encode, b0_encode, gradNonLin_encode, coilSense_encode, chemShift_encode, rot_encode, spatioTemporal_encode)
    %find which encode matrices have been selected:
    encode_matrices = {gradient_encode, offCtr_encode, b0_encode, gradNonLin_encode, coilSense_encode, chemShift_encode, rot_encode, spatioTemporal_encode};  % cell array of input matrices
    isEmptyFlags = cellfun(@isempty, encode_matrices);  % logical array: true if empty
    n = numel(encode_matrices);

    if isEmptyFlags(1) == 0
        ENCODE = gradient_encode; %gradient encode is always neccesary
    elseif isEmptyFlags(1) == 1
        error('You must have an encode matrix for gradient encoding. Please input a spiral trajectory.'); %edit error message later
    end 

    for i = 1:n
        if isEmptyFlags(n) == 0
            ENCODE = ENCODE.*encode_matrices(n);
        end
    end

end

