
function G = ft2(G)
G = fftshift(fft2(fftshift(G)));
end
