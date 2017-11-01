
function I2 = lensdistort(I,k)

[M N]=size(I); % Determine the size of the image to be distorted
center = [round(N/2) round(M/2)];

[xi,yi] = meshgrid(1:N,1:M); % Creates N x M (#pixels) x-y points

% Creates the mesh into a colum vector of coordiantes relative to the
% center

xt = xi(:) - center(1);
yt = yi(:) - center(2);

[theta,r] = cart2pol(xt,yt); % Converts the x-y coordinates to polar coordinates

R = sqrt(center(1)^2 + center(2)^2); % Calculate the maximum vector (image center to image corner)

r = r/R; % Normalize the polar coordinate r to range between 0 and 1 

s = r.*(1+k.*(r.^2)); % Aply the r-based transformation

s2 = s * R; % un-normalize s
  
if k <= 0 % Find a scaling parameter based on selected border type
    brcor = 1/(1 + k*(min(center)/R)^2);
    elseif k > 0
    brcor = r(1)/s(1);
end     
s2 = s2 * brcor;
[ut,vt] = pol2cart(theta,s2); % Convert back to cartesian coordinates
u = reshape(ut,size(xi)) + center(1);
v = reshape(vt,size(yi)) + center(2);
tmap_B = cat(3,u,v);
resamp = makeresampler('cubic','fill');
I2=tformarray(I,[],resamp,[2 1],[1 2],[],tmap_B,255);

end







