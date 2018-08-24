function [pattern_stack] = newEstimator(D1,PSF,shifts,approx_spacing,threshold)
frames = shifts^2;

N = length(D1(:,:,1));

%% Initialise all the values
%---------------------------------------------------------------------
PSF = PSF./max(PSF(:));
OTF_em = real(ift2(PSF));

N = length(PSF); 

D1 = D1 ./ max(D1(:));
wf1 = mean(D1,3);
wf1 = wf1./max(wf1(:));

%% Find maxima
%---------------------------------------------------------------------
% Fitting the grid to local mixima while odd has increased the program's
% ability to handle high background environments 
co_ords = zeros(2,1500,frames);
disp 'Finding maxima...'

for i = 1:frames
    b=find_maxima_gpu(D1(:,:,i),OTF_em,threshold);
    co_ords(:,1:length(b),i) = b;
end

%% Build maxima array
%---------------------------------------------------------------------
disp 'Building maxima_array...'

top = ceil((size(PSF,1))/2) - 15;
left = ceil((size(PSF,2))/2) - 15; 
bottom = top + 30;
right = left + 30;
small_PSF = PSF(top:bottom,left:right);
maxima_array = zeros(size(D1));
figure();
for i = 1:frames

     maxima_array(:,:,i) = zeros(N,N);
     for ii = 1:length(co_ords(:,:,i))
         if (co_ords(1,ii,i)+co_ords(2,ii,i))>=1
         top = co_ords(2,ii,i)-16;
         left = co_ords(1,ii,i)-16;
         bottom = top+30;
         right = left+30;
             if top>=1 && left>=1 && bottom<N && right<N
                  maxima_array(top:bottom,left:right,i) = maxima_array(top:bottom,left:right,i)+small_PSF;
             end
         end
     end      
end

%% Estimate spacings
%---------------------------------------------------------------------
disp 'Estimating spacings...'

P = 0;
av = [0 0; 0 0];
for k = 1:frames 
         
   if k ==1 || mod(k-1,shifts) == 0
   P = P+1;
   err = [0 0];
   for i = 1:10
      for ii = 1:10
         err(i,ii) = immse(circshift(maxima_array(:,:,k),[i+approx_spacing-5 ii+approx_spacing-5]),maxima_array(:,:,k));
      end
   end

   [error, shift] = min(err(:));
   [xshift, yshift] = ind2sub(size(err),shift);
   xshift = xshift+approx_spacing-4;
   yshift = yshift+approx_spacing-4;
   av(:,P) = [xshift yshift];
   disp(['Done frame: ' num2str(P) '/' num2str(shifts)])       
   
   end
end
av = sort(av);

av_x = ceil(mean(av(1,3:P-3))); % This excludes outliers and averages spacings
av_y = ceil(mean(av(2,3:P-3))); % This excludes outliers and averages spacings

disp (['x spacing = ' num2str(av_x)])
disp (['y spacing = ' num2str(av_y)])
pattern = generate_illum_patterns(N, PSF, av_x, av_y);

%% Estimate angle
%---------------------------------------------------------------------

disp 'Estimating angle...'
thetas = zeros(1,P);
P = 0;
pattern_stack = zeros(N,N,frames);
pattern = gpuArray((pattern));

% for k = 1:frames
% 
%    tx = av_x+2;
%    ty = av_y+2;
%    subx = 10;
%    suby = 10; 
%    if mod(k-1,shifts) == 0
%    full_pic = gpuArray(maxima_array(:,:,k));
%    P = P+1; 
%    a = 0;
%    err = gpuArray(zeros(tx,ty,6));   

%    for angle = -0.15:0.05:0.15
%        a = a+1;
%        for i = 1:tx
%           for ii = 1:ty
%                 temp = double(imrotate(single(pattern),angle,'crop','bicubic'));
%                 temp = circshift(temp,[i-subx ii-suby]);
%                 err(i,ii,a) = sqrt(sum(sum((temp(100:400,100:400)-full_pic(100:400,100:400)).^2)));
%           end
%        end
%    end
%    err = gather(err);
%    err(err<=0) = max(err(:));
%    [~, shift] = min(err(:));
%    [~, ~, inda] = ind2sub(size(err),shift);
%    thetas(P) = 0.05*(inda-1)-0.15;
%    disp (['Angle = ' num2str( thetas(P))])
%    disp(['Done frame: ' num2str(P) '/' num2str(shifts)])
%    end
% end
% 
% thetas = sort(thetas);
% av_theta = mean(thetas(3:P-3));
% disp (['Average angle = ' num2str(av_theta)])

%% Refind maxima
% Re-finding the spacings after rotation accounts for any increase/decrease
% due to the change in rotation
%---------------------------------------------------------------------
disp 'Building maxima_array...'
for i = 1:frames
%     D1(:,:,i) = imrotate(D1(:,:,i),-av_theta,'crop');
    b=find_maxima_gpu(D1(:,:,i),OTF_em,threshold);
    co_ords(:,1:length(b),i) = b;
end 

for i = 1:frames

      maxima_array(:,:,i) = zeros(N,N);
     for ii = 1:length(co_ords(:,:,i))
         if (co_ords(1,ii,i)+co_ords(2,ii,i))>=1
         top = co_ords(2,ii,i)-16;
         left = co_ords(1,ii,i)-16;
         bottom = top+30;
         right = left+30;
             if top>=1 && left>=1 && bottom<N && right<N
                  maxima_array(top:bottom,left:right,i) = maxima_array(top:bottom,left:right,i)+small_PSF;
             end
         end
     end
end

%% Re-estimate spacings
%---------------------------------------------------------------------
disp 'Estimating spacings...'
rav_x = round(av_x); rav_y = round(av_y);
P = 0;
av = [0 0; 0 0];
for k = 1:frames  
   if mod(k-1,shifts) == 0
   P = P+1;
   err = [0 0];
   for i = 1:6
      for ii = 1:6
         err(i,ii) = immse(circshift(maxima_array(:,:,k),[i+rav_x-3 ii+rav_y-3]),maxima_array(:,:,k));
      end
   end

   [~, shift] = min(err(:));
   [xshift, yshift] = ind2sub(size(err),shift);
   xshift = xshift+rav_x-3;
   yshift = yshift+rav_y-3;
   
   av(:,P) = [xshift yshift];
   disp(['Done frame: ' num2str(P) '/' num2str(shifts)])       
   
   end
end
av = sort(av);
av_x = (mean(av(1,3:P-3)));% This excludes outliers and averages spacings
av_y = (mean(av(2,3:P-3)));% This excludes outliers and averages spacings
disp (['final x-spacing = ' num2str(av_x)])
disp (['final y-spacing = ' num2str(av_y)])

pattern = generate_illum_patterns(N, PSF, av_x, av_y);

%% Determine shifts
%---------------------------------------------------------------------
disp 'Determining shifts...'
D1 = single(gpuArray(D1));
wf1 = single(gpuArray(mean(D1,3)));
step_size = av_x/sqrt(frames);
pattern_stack = single(gpuArray(repmat(pattern,[1 1 frames])));
clear maxima_array

% Generate shifted patterns

for i = 1:frames
    
x_shift = mod(i-1,sqrt(frames));
y_shift = floor(i/sqrt(frames));
x_shift = round(x_shift*step_size);
y_shift = round(y_shift*step_size);

pattern_stack(:,:,i) = (circshift(pattern_stack(:,:,i),[-x_shift -y_shift]));
end

%% Estimate shift
%--------------------------
x_counter = 0;
y_counter = 0;
err = single(gpuArray([0 0]));
for i = round(-av_x/2):round(av_x/2)
    x_counter = x_counter+1;
    y_counter = 0;
    for ii = round(-av_x/2):round(av_x/2)
        y_counter = y_counter+1;
        shifted = circshift(pattern_stack(:,:,1),[i ii]);
        shifted = shifted.*wf1;
        err(x_counter,y_counter) = sqrt(sum(sum((shifted-D1(:,:,1)).^2)));
    end
end

[~, shift] = min(err(:));
[xshift, yshift] = ind2sub(size(err),shift);
pattern_stack = (circshift(pattern_stack,[xshift-round(av_x/2) yshift-round(av_x/2)]));

%% Finalise shift
%--------------------------
tx = 10;
ty = 10;
err = single(gpuArray(zeros(4)));
for i = 1:tx
    for ii = 1:ty
        shifted = (circshift(pattern_stack,[ii-4 i-4]));
        shifted = pagefun(@times,shifted,wf1);
        shifted = pagefun(@minus,shifted,D1);
        shifted = pagefun(@times,shifted,shifted);
        err(i,ii) = sqrt(sum(shifted(:)));
    end
end

[~, shift] = min(err(:));
[xshift, yshift] = ind2sub(size(err),shift);
xshift = xshift-4;
yshift = yshift-4;
close all
pattern_stack=double(gather(pattern_stack));
D1 = double(gather(D1));
pattern_stack = circshift(pattern_stack,[xshift yshift]);
for p = 1:frames
    C = imfuse(pattern_stack(:,:,p),gather(D1(:,:,p)),'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    imshow(C); title 'initial overlay'; drawnow
end

disp(['Done frame: ' num2str(frames) '/' num2str(frames)]) 
close all
end
