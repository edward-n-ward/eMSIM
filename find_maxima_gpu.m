function [co_ords] = find_maxima_gpu(A,OTF_em,threshold)
A = gpuArray(A);
OTF_em = OTF_em./max(OTF_em(:));
dim = size(A);
OTF_em = imresize(OTF_em,[dim(1) dim(2)]);
A = real(ift2(ft2(A).*OTF_em));
A = A./max(A(:));
A(A<threshold) = 0;
C = imregionalmax(A);
C = imerode(C, ones(1));
[row col] = find(C);
dim = size(row);
co_ords = zeros(2, dim(1));
co_ords = gpuArray(co_ords);
co_ords(2,:) = row;
co_ords(1,:) = col;

co_ords(co_ords<10) = 0;
co_ords(co_ords>1014) = 0;

co_ords( ~any(co_ords,2), : ) = [];  %rows
co_ords( :, ~any(co_ords,1) ) = [];  %columns
co_ords = gather(co_ords);
end

