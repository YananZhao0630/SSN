function [G0,Gp,G_meta,G_radius] = Gabor_Bank_2D(N,radius,scales)
% GABOR_BANK_2D constructs a Gabor uniform covering frame in the Fourier
% domain 
%
% INPUT:
%  N                (pos int) frames are size NxN
%  radius           (pos int) F0 is supported in [-radius,radius]^2
%  scales           (pos int) number of Fourier scales
%
% OUTPUT:
%  G0               (array) lowest frequency frame
%  Gp               (array) higher frequency frames
%  G_meta           (array) frequency information of the frame
%  G_radius         (vector) radial support information
%--------------------------------------------------------------------------
% Weilin Li ~ June 2017

% initialize
G_num = ((2*scales+1)^2-1)/2;
Gp = zeros(N,N,G_num);
G_meta = zeros(2,G_num);
G_radius = zeros(1,scales);

% Fourier transform of f_0
G0 = Gabor_Eval_2D(N,radius,0,0); % create 2-D gabor atom

% Fourier transform of f_{m,n}
% frames are ordered clockwise with increasing radius

index = 0;
for s = 1:scales
    
    for n = -s:s
        index = index+1;
        Gp(:,:,index) = Gabor_Eval_2D(N,radius,s,n) + Gabor_Eval_2D(N,radius,-s,-n);
        G_meta(:,index) = [s;n];
    end
    
    for m = -(s-1):(s-1)
        index = index+1;
        Gp(:,:,index) = Gabor_Eval_2D(N,radius,-m,s) + Gabor_Eval_2D(N,radius,m,-s);
        G_meta(:,index) = [-m;s];
    end
    
    G_radius(s) = index;
end
    
