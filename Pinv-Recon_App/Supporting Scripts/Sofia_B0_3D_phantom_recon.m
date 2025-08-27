%% B0_phantom
%% Sofia Pearson
% 12/05/2024

%% import own data, or create "data" with phantom 
% data is created by multiplying ENCODE by phantom

%% make a switch case for 2D or 3D phantom
phantom_dimensions = '3D';


%% ============================================================================================================
%% FOR THE 3D CASE 
%  (eventually make this a switch case)

CondNumb=30; %what does this do??
ntx=64;%64,96,128

%% Generate k-space trajectory:
load('vd_spiral.mat'); %contains Kx0, Ky0, time0, (normal)
                      % Kx1, Ky1, time1 (undersampled)

%Check limits of k-space - ensure this is -0.5 to 0.5
% so then in gradient encode it gets x2pi and becomes -pi to pi?:
% Kx0 = (Kx0/(max(Kx0)))*0.5;
max(Kx0)
min(Kx0)

% Ky0 = (Ky0/(max(Ky0)))*0.5;
max(Ky0)
min(Ky0)

% plot k-space trajectory
% figure, plot(Kx0,Ky0), axis square, title('k-space trajectory'),

%% Generate 3D Shepp-Logan phantom:
addpath('C:\Users\sofia\OneDrive\Desktop\Pinv-Recon_Dev\misc\public\dependencies\')
Phantom3D = phantom3d('Modified Shepp-Logan',ntx); 

nslices = 1; %3 % Choose how many head slices you want to see

%remove all slices we are not imaging:
% such that dimensions of 3d phantom go from 
% ntx*ntx*ntx to ntx*ntx*nslices
chosen_slices = zeros(nslices, 1);
for n = 1:nslices
    slice = round((n/(nslices+1))*ntx);
    chosen_slices(n,1) = slice;
end

Phantom3D = Phantom3D(:,:,chosen_slices);
figure, mat2montage(Phantom3D), title('Chosen Slices');


%% Generate the B0 maps: 
x = linspace(-1,1,ntx);
y = x;
[X, Y] = ndgrid(x, y);

% figure, imagesc(X),title('X')
% figure, imagesc(Y),title('Y')

B0Map=125*Y.^2-30; % why are these scale factors chosen?
% field map range across 125 Hz total variation with a -30 Hz global shift at the center
% figure, imagesc(B0Map), title('B0 map')

%% Gradient Encoding
Kx=Kx0; Ky=Ky0; time=time0; % Kx and ky are just defining your k-space spiral
Gradient_Encode = single(exp(1i*2*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).')));
%^^ this line is basically just encoding the k-space trajectory 

%if k space is from -1 to 1
%Gradient_Encode = single(exp(1i*pi*(Kx(:)*X(:).'+Ky(:)*Y(:).'))); 

%% Off Resonance Encoding
OffRes_Encode = single(exp(1i*time(:)*(B0Map(:)).')); % if units is in Herts

%OffRes_Encode = single(exp(1i*2*pi*time(:)*(B0Map(:)).')); % if units are
%in radians per second

%% Get overall Encode Matrix
Encode = Gradient_Encode .* OffRes_Encode;

% now get all the plots and do pinv recon for each plot:
for n = 1:nslices

    %==============================================
    fprintf('>> Computing reconstructions for Slice %d of %d ...', n, nslices);
    %==============================================
        % disp, step 1 of n
    imgPhantom = squeeze(Phantom3D(:,:,n));
    
    %% Get Gradient Encode Only Reconstruction
    tic, [U1,S1,V1]=svd((Gradient_Encode),'econ'); toc, 
    diagS1=diag(S1); imax=find(diag(S1)>max(diag(S1))/CondNumb,1,'last'); invS1=1./diag(S1); invS1(imax+1:end)=0; invS1=diag(invS1);
    
    tic, Gradient_Recon =V1*invS1*U1'; toc, 
    
    Img_Gradient_only = reshape(Gradient_Recon*(Encode*imgPhantom(:)), [ntx, ntx]);
    %^ Encode*imgPhantom2(:) created synthetic "data"

    %% Get Grad + OffRes Encoding Reconstruction
    tic, [U,S,V]=svd((Encode),'econ'); toc, 
    diagS=diag(S); imax=find(diag(S)>max(diag(S))/CondNumb,1,'last'); invS=1./diag(S); invS(imax+1:end)=0; invS=diag(invS);
    
    tic, Recon =V*invS*U'; toc, 
    
    IMG = reshape(Recon*(Encode*imgPhantom(:)), [ntx, ntx]);

    %% Plot Results
    figure, tiledlayout(1,2)

    nexttile,
    imagesc(abs(Img_Gradient_only)), title('Pinv'), axis square, 
    
    nexttile,
    imagesc(abs(IMG)), title('Pinv + OffRes'), axis square, 

end

