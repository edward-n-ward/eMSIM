function [pattern]=generate_illum_patterns2(N,PSF_ex1,S)

K = S^2;
%-----------------------------------------------------
% Generate PSF patterns P at the sample plane 
%-----------------------------------------------------

top = ceil((size(PSF_ex1,1))/2) - 10;
left = ceil((size(PSF_ex1,2))/2) - 10; 
bottom = top + 20;
right = left + 20;
small_PSFex1 = zeros(20,20);
PSF = PSF_ex1(top:bottom,left:right);
PSF = PSF./max(PSF(:));

%-----------------------------------------------------
% Build pattern 
%-----------------------------------------------------

PSF_sizex = size(PSF,1);
PSF_sizey = size(PSF,2);
pImageData = zeros(N+3*PSF_sizex,N+3*PSF_sizey);
shifts = floor(N/S);
%%
% iterate over all of the spots

for ii = 1:S:N+2*PSF_sizex
    for i = 1:S:N+2*PSF_sizex

        top = round(ii);
        bottom = top+PSF_sizey-1;
        left = round(i);
        right = left+PSF_sizex-1;

        pImageData(top:bottom,left:right) = pImageData(top:bottom,left:right) + PSF;
        
    end
end
%figure(); imagesc(pImageData);

pattern = zeros(N,N,K);
k = 0; 
for l = 1:S
    pImageData = circshift(pImageData,[0 -1]);
    temp = pImageData;
    for m = 1:S
        k = k+1;        
        temp = circshift(temp,[-1 0]);  
        p = floor(PSF_sizex/2);
        n = floor(PSF_sizey/2);
        pattern(:,:,k) = temp(p:N-1+p,n:N-1+n);
    end
end

end
