% SCAT_COMPARE computes the features produced by the windowed
% scattering transform and the Fourier scattering transform.
% This script requires ScatNet, which can be downloaded here
% http://www.di.ens.fr/data/
% This script reproduces the numerical experiment in the paper,
% "Analysis of time-frequency scattering transforms",
% by Wojciech Czaja and Weilin Li
%--------------------------------------------------------------------------
% Weilin Li ~ June 2017
clc;clear;
%% Initialization 
path_to_scatnet = fileparts(mfilename('fullpath'));

addpath(fullfile(path_to_scatnet, 'Scattering'));
addpath(fullfile(path_to_scatnet, 'Utils'));
addpath(fullfile(path_to_scatnet, 'Utilities'));

% sample image
f = double(imread('lena_gray_512.tif'));
f = f/norm(f,'fro');

% threshold cutoff
cutoff = norm(f,'fro')*.005;% A和A‘的积的对角线和的平方根，即sqrt(sum(diag(A'*A)))

%% Fourier Scattering Transform

% set parameters
radius = 24;% atom size
width = 2;  % scale distribution
depth = 2;  % layer
dec = 'yes'; %option for path decreasing coefficients

% compute the FST coefficients
[FST,FST_meta] = Fourier_Scattering_2D(f,radius,width,depth,dec);

% threshold coefficients tensor
% discard the coefficients which is smaller
[FST_thres, FST_meta_thres] = threshold_coeff(FST,FST_meta,cutoff);

%% add energy
size(FST{1,2})
X_en = sum(sum(f.^(2)));
% m = 0 时的能量
X_S_0 = sum(sum((FST{1}).^(2)));
ratio_0 = sum((X_S_0./X_en)*100)



%% Wavelet Scattering Transform 

% set parameters
filt_opt.J = 2;
filt_opt.L = 10;
scat_opt.M = 2;

% compute coefficients
[Wop, filters] = wavelet_factory_2d(size(f), filt_opt, scat_opt);
WST = scat(f,Wop);
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
end

%% create aggregate image 

N = size(WST_thres,1);
ds = size(FST_thres,1)/size(WST_thres,1); % resize factor
K = 10;
FST_num = size(FST_thres,3);
WST_num = size(WST_thres,3);
FST_im = zeros(N*ceil(FST_num/K),K*N);
WST_im = zeros(N*ceil(FST_num/K),K*N);

for j = 1:ceil(FST_num/K)
    for k = 1:K
        m = (j-1)*K+k;
        if m <= FST_num
            temp = FST_thres(:,:,m);
            temp = temp(1:ds:end,1:ds:end);
            FST_im(1+(j-1)*N:j*N,1+(k-1)*N:k*N) = temp/max(max(temp));
        else
            FST_im(1+(j-1)*N:j*N,1+(k-1)*N:k*N) = ones(N,N);
        end
    end
end

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
end

figure; imshow([FST_im;WST_im])
size(WST{1,2})
X_en = sum(sum(f.^(2)));
% m = 0 时的能量
X_S_0 = sum(sum((WST{1,2}.signal{1,1}).^(2)));
ratio_0 = sum((X_S_0./X_en)*100)
