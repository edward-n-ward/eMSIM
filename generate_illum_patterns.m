function [pImageData, PSF_sizex] = generate_illum_patterns(N,PSF, xshift, yshift)

top = ceil((size(PSF,1))/2)-ceil((size(PSF,1))/16);
left = ceil((size(PSF,2))/2)-ceil((size(PSF,2))/16); 
bottom = top+ceil((size(PSF,1))/8);
right = left+ceil((size(PSF,1))/8);
PSF = PSF(top:bottom,left:right);

%%
%calculate parameters
PSF_sizex = size(PSF,1);
PSF_sizey = size(PSF,2);
pImageData = zeros(N+3*PSF_sizex,N+3*PSF_sizey);

%%
% iterate over all of the spots
for i = 1:xshift:N+2*PSF_sizex
    for ii = 1:yshift:N+2*PSF_sizey
    top = ii;
    bottom = top+PSF_sizey-1;
    left = i;
    right = left+PSF_sizex-1;
        if size(pImageData(top:bottom,left:right)) == size(PSF)
        pImageData(top:bottom,left:right) = pImageData(top:bottom,left:right) + PSF;
        end
    end
end

end