tic
clearvars
close all
%% Input parameters
%-----------------------------------------------------
Read_out_noise_on = 0;
S = 15; % pitch of the square pattern
K = S^2;
piFP_iter = 1; % number of piFP iterations
jRL_iter = 1; % number of jRL iterations
sigma = 1e-11; % variance of gaussian noise
threshold = 0.4; % threshold for peak finding
subtract_const = 0.5; %threshold for difference subtraction

%% Create save location
%-----------------------------------------------------
%the data files can get very large so to save RAM we put
%them in a temporary file on the desktop
d_top = winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop');
location = strcat(d_top,'\temp image bin','.mat');
m=matfile(location,'writable',true);

%% Generate sample
%-----------------------------------------------------
% Reads the sample from file
Sample = single(imread('seimens.tif', 'tif'));
Sample = mean(Sample,3);
Sample = Sample./max(Sample(:));
N = size(Sample,1); % number of pixels

%% Read optical functions of the simulated system
%-----------------------------------------------------
disp('Running...');
PSF_ex1 = (single(imread('PSF_ex.tif', 'tif')))/255;
PSF_ex2 = (single(imread('spiral phase 2.tif', 'tif')))/255;
PSF_em = (single(imread('PSF_em.tif', 'tif')))/255;
OTF_em = single(real(ift2(PSF_em)));

%% Generate illumination patterns
%-----------------------------------------------------
disp('Generating illumination patterns...');
top = ceil((size(PSF_ex1,1))/2) - 10;
left = ceil((size(PSF_ex1,2))/2) - 10; 
bottom = top + 20;
right = left + 20;
small_PSFex1 = PSF_ex1(top:bottom,left:right);
small_PSFex2 = PSF_ex2(top:bottom,left:right);

P1 = generate_illum_patterns2(N, small_PSFex1, S);
P1 = single(P1);
P2 = generate_illum_patterns2(N, small_PSFex2, S);
P2 = single(P2);
P1 = P1./max(P1(:));
P2 = P2./max(P2(:));
m.P2 = P2;
clear P2
clear  top bottom left right;

%% Generate 1st Raw data 
%-----------------------------------------------------
D1 = zeros(N,N,K,'single');
Dgpu = gpuArray(D1);
P1 = gpuArray(P1);
disp('Generating raw data 1...')

for i = 1:K
    Dgpu(:,:,i) = ift2(ft2(P1(:,:,i).*Sample).*OTF_em);
    Dgpu(:,:,i) = abs(Dgpu(:,:,i));
end
P1 = gather(P1);
m.P1 = P1;
clear P1
if Read_out_noise_on 
    disp('Adding Gaussian noise...');
    for i = 1:K
        Dgpu(:,:,i) = imnoise(Dgpu(:,:,i),'gaussian', 0,sigma);
    end
end
D1 = gather(Dgpu);
D1 = D1./max(D1(:));
wf_down = mean(D1,3);
% figure(); imagesc(wf_down(225:275,325:375)); axis square; axis off; colormap gray; title('Widefield 2 zoom');
m.D1 = D1;
clear D1 
%% Generate 2nd Raw data 
%-----------------------------------------------------
disp 'Generating raw data 2...';
P2 = m.P2;
P2 = gpuArray(P2);
for i = 1:K
    Dgpu(:,:,i) = ift2(ft2(P2(:,:,i).*Sample).*OTF_em);
    Dgpu(:,:,i) = abs(Dgpu(:,:,i));
end
clear sample
clear P2
if Read_out_noise_on 
    disp 'Adding Gaussian noise 2...';
    for i = 1:K
        Dgpu(:,:,i) = imnoise(Dgpu(:,:,i),'gaussian', 0,sigma);
    end
end

D2 = gather(Dgpu);
clear Dgpu
D2 = D2./max(D2(:));
m.D2 = D2;
clear D2 P2

%% Estimate patterns
%---------------------------------------------------------------------
 disp 'Estimating patterns...'
%  un-comment this code if you wish to use the estimator function
%  D1 = m.D1;
%  D1 = double(D1);
%  P = newEstimator(D1,PSF_ex1,S,S,threshold);
%  D1 = single(D1);
%  m.D1 = D1;
%  clear D1
 co_ords = zeros(2,1500,K);
 P = m.P1;
 for i = 1:K
 b=find_maxima_gpu(P(:,:,i),OTF_em,threshold);
 co_ords(:,1:length(b),i) = b;
 end
 clear b
 P = single(P);
 m.P = P;
 clear P
 
%% Build doughnut array
 %---------------------------------------------------------------------
 disp 'Building doughnut array...'
 P2 = zeros(N,N,K);
 for i = 1:K
     for ii = 1:length(co_ords(:,:,i))
     if (co_ords(1,ii,i)+co_ords(2,ii,i))>=1
     top = co_ords(2,ii,i)-11;
     left = co_ords(1,ii,i)-11;
     bottom = top+20;
     right = left+20;
         if top>=1 && left>=1 && bottom<N && right<N
             P2(top:bottom,left:right,i) = P2(top:bottom,left:right,i)+small_PSFex2;
         end
     end
     end
 end

m.P2 = P2;
clear P2

%% Generate IWS pattern
%---------------------------------------------------------------------
disp 'Generating IWS pattern...'
psf_IWS = (small_PSFex1-subtract_const*small_PSFex2);
psf_IWS(psf_IWS<0)=0;
P_IWS = single(zeros(N,N,K));
 for i = 1:K
     for ii = 1:length(co_ords(:,:,i))
     if (co_ords(1,ii,i)+co_ords(2,ii,i))>=1
     top = co_ords(2,ii,i)-11;
     left = co_ords(1,ii,i)-11;
     bottom = top+20;
     right = left+20;
         if top>=1 && left>=1 && bottom<N && right<N
            P_IWS(top:bottom,left:right,i) = P_IWS(top:bottom,left:right,i)+psf_IWS;
         end
     end
     end
 end
 m.P_IWS = P_IWS;
 clear P_IWS 
 
%% Generate IWS data
%---------------------------------------------------------------------
disp 'Generating IWS data...'
D1 = m.D1;
D2 = m.D2;
 D_IWS = single(zeros(size(D1)));
 for i = 1:K
     temp = D1(:,:,i)-subtract_const*D2(:,:,i);
     temp(temp<0) = 0;

     D_IWS(:,:,i) = temp;
 end
 D_IWS = D_IWS./max(D_IWS(:));
 clear D1 D2
 m.D_IWS = D_IWS;
  clear D_IWS
 
 %% Upsample
%---------------------------------------------------------------------
% This upsample section is not needed given sufficient resolution 
% of the PSF image file
disp 'Resizing...'
PSF_ex1 = interp2(PSF_ex1);
Sample=interp2(Sample);
%---------------------------------------------------------------------
D1 = m.D1;
D1_gpu = gpuArray(D1); 
clear D1
for i = 1:K
temp = interp2(D1_gpu(:,:,i));
D3(:,:,i) = gather(temp);
end
clear D1_gpu 
wf1 = mean(D3,3);
m.wf1 = wf1;
clear wf1
m.D1 = D3;

%---------------------------------------------------------------------
D2 = m.D_IWS;
D2_gpu = gpuArray(D2); 
clear D2
for i = 1:K
temp = interp2(D2_gpu(:,:,i));
D3(:,:,i) = gather(temp);
end
clear D2_gpu 
m.D_IWS = D3;
FED_image = mean(D3,3);
%  figure();imagesc(FED_image(400:600,400:600));colormap gray;axis square; title 'FED image zoom'
%  figure();imagesc(FED_image);colormap gray;axis square; title 'FED image'
 temp = double(FED_image);
 temp = temp./max(temp(:));
 imwrite(temp,'FED image.tif')
 m.FED = FED_image;
clear D3 FED_image
%---------------------------------------------------------------------
P2 = m.P1;
P_gpu = gpuArray(P2); 
clear P2
for i = 1:K
temp = interp2(P_gpu(:,:,i));
P2(:,:,i) = gather(temp);
end
clear P_gpu 
m.P1 = P2;
clear P2
%---------------------------------------------------------------------
P = m.P_IWS;
P_gpu = gpuArray(P); 
clear P
for i = 1:K
temp = interp2(P_gpu(:,:,i));
P(:,:,i) = gather(temp);
end
clear P_gpu temp 
m.P_IWS = P;
clear P

OTF_em = interp2(OTF_em);


%% MSIM reconstruction piFP
%---------------------------------------------------------------------
% We have not yet found a way to implement the PiFP on the GPU as
% the necessary stacks take up too much memory
D_IWS = m.D1;
wf1 = mean(D_IWS,3);
P_IWS = m.P1;
Ig = wf1;
if piFP_iter > 0
disp('Reconstruction using pattern-illuminated Fourier Ptychography...');
FD = zeros(size(D_IWS),'single');
for ii = 1:K
    FD(:,:,ii) = ft2(D_IWS(:,:,ii));
end
clear D_IWS

p = randperm(K);

    for n=1:piFP_iter 
         disp(['piFP iteration ' num2str(n)]); 
         for ii = 1:K
            frame = p(ii);
            P_ill = P_IWS(:,:,frame);
            It = Ig.*P_ill;
            ft_It = ft2(It);
            ft_It = ft_It + conj(OTF_em).*(FD(:,:,frame) - OTF_em.*ft_It); 
            It_upd = abs(ift2(ft_It));
            Ig_upd = Ig + P_ill./(max(P_ill(:))^2).*(It_upd - It);  % object update      
            Ig = Ig_upd;
         end          
    end
clear FD Ig_upd It It_upd p P_ill
end
%% MSIM Reconstruction jRL
%---------------------------------------------------------------------

disp('Reconstruction using pattern-illuminated jRL...');
D_IWS = m.D1; 
D_IWS = gpuArray(D_IWS);
P_IWS = gpuArray(P_IWS);
for n=1:jRL_iter
    disp(['jRL iteration ' num2str(n)]);
    intermediate = zeros(size(wf1));     

     for ii = 1:K
         
        temp = Ig.*P_IWS(:,:,ii);
        temp = ft2(temp);
        temp = abs(ift2(temp.*OTF_em)); % modeled images                  
        temp = (D_IWS(:,:,ii))./(temp);                    
        temp = ft2(temp);
        temp = abs(ift2(temp.*OTF_em));
        intermediate = intermediate + temp.*P_IWS(:,:,ii);

     end

     Ig = Ig.*intermediate; 
     
end     
clear temp P_IWS D_IWS
% figure();imagesc(wf1);colormap gray;axis square; title 'Widefield';
% figure();imagesc(wf1(400:600,400:600));colormap gray;axis square; title 'Widefield zoom'
% figure();imagesc(Ig);colormap gray;axis square; title 'MSIM reconstruction';
% figure();imagesc(Ig(400:600,400:600));colormap gray;axis square; title 'MSIM reconstruction zoom';

m.MSIM = Ig;
%% eMSIM Reconstruct piFP
%---------------------------------------------------------------------
D_IWS = m.D_IWS;
wf2 = mean(D_IWS,3);
P_IWS = m.P_IWS;
if piFP_iter > 0
disp('Reconstruction using pattern-illuminated Fourier Ptychography...');
FD = zeros(size(D_IWS),'single');
for ii = 1:K
    FD(:,:,ii) = ft2(D_IWS(:,:,ii));
end
clear D_IWS

p = randperm(K);
Ig = wf2;

    for n=1:piFP_iter 
         disp(['piFP iteration ' num2str(n)]); 
         for ii = 1:K
            frame = p(ii);
            P_ill = P_IWS(:,:,frame);
            It = Ig.*P_ill;
            ft_It = ft2(It);
            ft_It = ft_It + conj(OTF_em).*(FD(:,:,frame) - OTF_em.*ft_It); 
            It_upd = abs(ift2(ft_It));
            Ig_upd = Ig + P_ill./(max(P_ill(:))^2).*(It_upd - It);  % object update      
            Ig = Ig_upd;
         end          
    end
clear FD Ig_upd It It_upd p P_ill
end
%% Reconstruct jRL
%---------------------------------------------------------------------

disp('Reconstruction using pattern-illuminated jRL...');
D_IWS = m.D_IWS; 
D_IWS = gpuArray(D_IWS);
P_IWS = gpuArray(P_IWS);

for n=1:jRL_iter
    disp(['jRL iteration ' num2str(n)]);
    intermediate = zeros(size(wf1));     

     for ii = 1:K
         
        temp = Ig.*P_IWS(:,:,ii);
        temp = ft2(temp);
        temp = abs(ift2(temp.*OTF_em)); % modeled images                  
        temp = (D_IWS(:,:,ii))./(temp);                    
        temp = ft2(temp);
        temp = abs(ift2(temp.*OTF_em));
        intermediate = intermediate + temp.*P_IWS(:,:,ii);

     end

     Ig = Ig.*intermediate; 
     
end

clear temp P_IWS D_IWS

%% Display Results
%---------------------------------------------------------------------
Ig = double(gather(Ig));
Ig = Ig./max(Ig(:));
% figure();imagesc(Ig);colormap gray;axis square; title 'eMSIM reconstruction zoom'
MSIM = m.MSIM;
FED = m.FED;
