function [pattern_stack] = rebuilt_estimator(D1,PSF,shifts,approx_spacing,threshold)

%% Initialise all the values
%---------------------------------------------------------------------
PSF = PSF./max(PSF(:));
OTF_em = abs(ft2(PSF));
frames = shifts^2; % How many images in the stack to find the pattern for
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
 %Read the appropriate image from the stack
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
   yshift = yshift+approx_spacing-5;
   av(:,P) = [xshift yshift];
   disp(['Done frame: ' num2str(P) '/' num2str(shifts)])       
   
   end
end
av = sort(av);

av_x = ceil(mean(av(1,3:P-3))); % This excludes outliers and averages spacings
av_y = ceil(mean(av(2,3:P-3))); % This excludes outliers and averages spacings

disp (['x spacing = ' num2str(av_x)])
disp (['y spacing = ' num2str(av_y)])
pattern = generate_illum_patterns_small(N, PSF, av_x, av_y);

%% Estimate angle
%---------------------------------------------------------------------
% Put the lens distortion correction in here
disp 'Estimating angle...'
thetas = zeros(1,P);
P = 0;
pattern_stack = zeros(N,N,frames);
pattern = gpuArray((pattern));

for k = 1:frames

   tx = av_x+2;
   ty = av_y+2;
   subx = 10;
   suby = 10; 
   if mod(k-1,shifts) == 0
   full_pic = gpuArray(maxima_array(:,:,k));
   P = P+1; 
   a = 0;
   err = gpuArray(zeros(tx,ty,6));   

   for angle = -0.3:0.1:0.3
       a = a+1;
       for i = 1:tx
          for ii = 1:ty
                temp = double(imrotate(single(pattern),angle,'crop'));
                temp = circshift(temp,[i-subx ii-suby]);
                err(i,ii,a) = sqrt(sum(sum((temp-full_pic).^2)));
          end
       end
   end
   err = gather(err);
   err(err<=0) = max(err(:));
   [~, shift] = min(err(:));
   [~, ~, inda] = ind2sub(size(err),shift);
   thetas(P) = 0.1*(inda-1)-0.2;
   disp (['Angle = ' num2str( thetas(P))])
   disp(['Done frame: ' num2str(P) '/' num2str(shifts)])
   end
end

thetas = sort(thetas);
av_theta = mean(thetas(3:P-3));
disp (['Average angle = ' num2str(av_theta)])

%% Refind maxima
% Re-finding the spacings after rotation accounts for any increase/decrease
% due to the change in rotation
%---------------------------------------------------------------------
disp 'Building maxima_array...'
for i = 1:frames
    D1(:,:,i) = imrotate(D1(:,:,i),-av_theta,'crop');
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

P = 0;
av = [0 0; 0 0];
for k = 1:frames  
   if mod(k-1,shifts) == 0
   P = P+1;
   err = [0 0];
   for i = 1:6
      for ii = 1:6
         err(i,ii) = immse(circshift(maxima_array(:,:,k),[i+av_x-3 ii+av_y-3]),maxima_array(:,:,k));
      end
   end

   [~, shift] = min(err(:));
   [xshift, yshift] = ind2sub(size(err),shift);
   xshift = xshift+av_x-2;
   yshift = yshift+av_y-3;
   av(:,P) = [xshift yshift];
   disp(['Done frame: ' num2str(P) '/' num2str(shifts)])       
   
   end
end
av = sort(av);
av_x = ceil(mean(av(1,3:P-3)));% This excludes outliers and averages spacings
av_y = ceil(mean(av(2,3:P-3)));% This excludes outliers and averages spacings

disp (['final x-spacing = ' num2str(av_x)])
disp (['final y-spacing = ' num2str(av_y)])
pattern = generate_illum_patterns_small(N, PSF, av_x, av_y);

%% Determine shifts
%---------------------------------------------------------------------
disp 'Determining shifts...'

tx = av_x+5;
ty = av_y+5;
subx=round(av_x/2); 
suby=round(av_y/2);
for k = 1:frames
    % Can we constrain vertical shifts after the first frame?
    if  mod(k-1,shifts) == 0 
        err = gpuArray([0 0]);
        pattern2 = gpuArray(pattern);
        temp_gpu = gpuArray(maxima_array(:,:,k));
        for i = 1:(tx)
            for ii = 1:(ty)
             temp = circshift(pattern2,[i-subx ii-suby]);
             err(i,ii) = sqrt(sum(sum((temp_gpu-temp).^2)));
            end
        end
        [~, shift] = min(err(:));
        [xshift, yshift] = ind2sub(size(err),shift);
        xshift = xshift-subx;
        yshift = yshift-suby;
        disp(['Done frame: ' num2str(k-1) '/' num2str(frames)]) 
         disp(['x-shift: ' num2str(xshift)])
         disp(['y-shift: ' num2str(yshift)])
        pattern2 = gather(pattern2);
    else
        pattern2 = pattern_stack(:,:,k-1);
        err = [0 0];
        for i = 1:12
            for ii = 1:12
             temp = circshift(pattern2,[i-6 ii-6]);
             err(i,ii) = immse(maxima_array(:,:,k),temp);
            end
        end
        [~, shift] = min(err(:));
        [xshift, yshift] = ind2sub(size(err),shift);
        xshift = xshift -6;
        yshift = yshift -6;
        disp(['x-shift: ' num2str(xshift)])
        disp(['y-shift: ' num2str(yshift)])
    end
    pattern_stack(:,:,k) = circshift(pattern2,[xshift yshift]);
end
disp(['Done frame: ' num2str(frames) '/' num2str(frames)]) 
end
