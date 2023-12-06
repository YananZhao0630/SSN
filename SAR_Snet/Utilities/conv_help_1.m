function U = conv_help_1(f,Gp)
% CONV_HELP_1 computes absolutle value of the convolution of f and Gp
%
% INPUT:
%  f                (array) function in spatial domain
%  Gp               (tensor) frames in Fourier domain
%
% OUTPUT:
%  U                (tensor) the absolute value of the convolution
%--------------------------------------------------------------------------
% Weilin Li ~ June 2017

F = fftshift(fft2(f));
U = zeros(size(Gp));
for j = 1:size(Gp,3)
    U(:,:,j) = abs(ifft2(ifftshift(Gp(:,:,j).*F)));
end

