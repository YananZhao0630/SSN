function S = conv_help_2(G0,U)
% CONV_HELP_2 computes the real part of the convolution of G0 and U
%
% INPUT:
%  G0               (array) frame in Fourier domain
%  U                (tensor) functions in spatial domain
%
% OUTPUT:
%  S                (tensor) of convolutions
%--------------------------------------------------------------------------
% Weilin Li ~ June 2017

S = zeros(size(U));
for j = 1:size(U,3)
    S(:,:,j) = real(ifft2(ifftshift(fftshift(fft2(U(:,:,j))).*G0)));
end

