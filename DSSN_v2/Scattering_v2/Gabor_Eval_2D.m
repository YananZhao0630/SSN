function G = Gabor_Eval_2D(N,A,m,n)
% GABOR_EVAL evaluates the Fourier transform of g_{m,n}

% INPUT:
%  N                (pos int) each frame is a matrix of size NxN
%  A                (pos int) supported in a square of length 2A
%  m                (int) frame index
%  n                (int) frame index
%
% OUTPUT:
%  G                (matrix) of function values at the grid points (x,y)
%--------------------------------------------------------------------------
% Weilin Li ~ June 2017

[x,y] = meshgrid(-N/2+1/2:N/2-1/2);
cutoff = (abs(x-m*A) < A) .* (abs(y-n*A) < A);%%去除边界以外的
G = cutoff .* cos(pi*(x-m*A)/(2*A)) .* cos(pi*(y-n*A)/(2*A));
