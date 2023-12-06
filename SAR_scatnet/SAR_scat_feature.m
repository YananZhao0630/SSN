% SCAT_COMPARE computes the features produced by the windowed
% scattering transform and the Fourier scattering transform.
% This script requires ScatNet, which can be downloaded here
% http://www.di.ens.fr/data/
% This script reproduces the numerical experiment in the paper,
% "Analysis of time-frequency scattering transforms",
% by Wojciech Czaja and Weilin Li
%--------------------------------------------------------------------------
% Weilin Li ~ June 2017

%% Initialization 
%clc;clear;
% sample image
clc;clear;
addpath('./Examples');
addpath('./Images');
addpath('./Scattering');
addpath('./Utilities');
addpath('./Utils');
addpath('./cwnn');


f1 = double(imread('Sulzberger1_1.bmp'))/256;
f2 = double(imread('Sulzberger1_2.bmp'))/256;

f_gray1 = rgb2gray(f1);
f_gray2 = rgb2gray(f2);
%imagesc(f_gray)
f_gray1 = f_gray1/norm(f_gray1,'fro');
f_gray2 = f_gray2/norm(f_gray2,'fro');
f_gray = f_gray1./f_gray2;

f_gray = f_gray/norm(f_gray,'fro');
% threshold cutoff
cutoff = norm(f_gray,'fro')*.005;

%% Fourier Scattering Transform


% set parameters
radius = 24;
width = 4;
depth = 2;
dec = 'yes';



%% Windowed Scattering Transform 

% set parameters
filt_opt.J = 4;
filt_opt.L = 10;
scat_opt.M = 1;
scat_opt.oversampling = 4;



% compute coefficients
[Wop, filters] = wavelet_factory_2d(size(f_gray), filt_opt, scat_opt);


[WST,newU,newV] = scat(f_gray,Wop);
[WST_mat,WST_meta] = format_scat(WST);
WST_mat = double(WST_mat);

% find largest coefficients
WST_j_thres = [];
WST_theta_thres = [];
WST_thres = [];
for k = 1:size(WST_mat,1)
    if norm(squeeze(WST_mat(k,:,:)),'fro') >= cutoff
        WST_thres = cat(3,WST_thres,squeeze(WST_mat(k,:,:)));
        WST_j_thres = [WST_j_thres, WST_meta.j(:,k)];
        WST_theta_thres = [WST_theta_thres, WST_meta.theta(:,k)];
    end
    k
end

%% create aggregae image 

N = size(WST_thres,1);
K = 14;
WST_num = size(WST_thres,3);
WST_im = zeros(N,K*N);



for j = 1:ceil(WST_num/K)
    for k = 1:K
        m = (j-1)*K+k;
        if m <= WST_num
            temp = WST_thres(:,:,m);
            WST_im(1+(j-1)*N:j*N,1+(k-1)*N:k*N) = temp/max(max(temp));
        else
            WST_im(1+(j-1)*N:j*N,1+(k-1)*N:k*N) = ones(N,N);
        end
    end
    j
end

figure; imshow([WST_im])
%imshow([abs(FST_im(1:64,65:128))-FST_im(1:64,65:128)])

