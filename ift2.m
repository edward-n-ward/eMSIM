
function G = ift2(G)
G = ifftshift(ifft2(ifftshift(G)));
end