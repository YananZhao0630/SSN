function xdown = mydownsample(x,r)
% downsamples x by a factor of r
% if x is a matrix, downsamples each column by a factor of r

xdown = [];
if ismatrix(x)
    for jj = 1:size(x,2)
        xdown(:,jj) = x(1:r:end, jj);
    end
    
elseif isvector(x)
    xdown = x(1:r:end); 
    
else 
    error('datatype not supported')
end
