function [pattern, k] = optimiseskew(pattern,image)
f = @(k,pattern,image)immse(lensdistort(pattern,k),image);
fun = @(k)f(k,pattern,image);
estimate = 0;
k = fminsearch(fun,estimate);
pattern = lensdistort(pattern,k);
end
