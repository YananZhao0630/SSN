function [FST,FST_meta] = Fourier_Scattering_2D(f,radius,scales,depth,dec)
% FOURIER_SCATTERING_2D computes the truncated Fourier scattering
% transform, see GABOR_BANK_2D
%
% INPUT:
%  f                (array) square image
%  radius           (pos int) frame parameter
%  width            (pos int) approximately the number of childern per node
%  depth            (pos int) number of layers
%  dec              (string) option for path decreasing coefficients
%
% OUTPUT:
%  FST              (cell) Fourier scattering coefficients
%  FST_meta         (cell) frequency information of the coeffcients
%--------------------------------------------------------------------------
% Weilin Li ~ June 2017

% initialize
U = cell(1,depth+1);
FST = cell(1,depth+1);
FST_meta = cell(1,depth+1);

% construct Gabor frame
[G0,Gp,G_meta,G_radius] = Gabor_Bank_2D(size(f,1),radius,scales);

% zero order coefficients
U{1} = f;
FST{1} = real(ifft2(ifftshift(fftshift(fft2(f)).*G0)));%F-1（F(f)*G）

% for each layer
for n = 1:depth    
    
    % for each function in the current layer
    for j = 1:size(U{n},3) % number of slices at each layer
        
        % fix a coefficient
        f_cur = U{n}(:,:,j);
        
        % determine the Fourier support of the coefficient
        if n == 1
            f_meta = [];
            f_freq = scales;
        elseif strcmp(dec,'yes') 
            f_meta = FST_meta{n}(:,j);
            f_freq = max(abs(f_meta(2*(n-1)-1:2*(n-1))))-1;
        else
            f_meta = FST_meta{n}(:,j);
            f_freq = scales;
        end
        
        % check if the frequency is positive
        % plus and minus is common in one gabor
        if f_freq > 0 
            
            % absolute value of the convolution of f_cur with the frames
            G_num = G_radius(f_freq);
            U_next = conv_help_1(f_cur,Gp(:,:,1:G_num));
            % F-1（F(f_cur)*G）
            S_next = conv_help_2(G0,U_next);
            % real(F-1(F(U_next)*G0))
            
            % update Scat, Ucat, and Meta
            U{n+1} = cat(3,U{n+1},U_next);% sum the input of the next layer
            % C = cat(dim, A, B) 按dim来联结A和B两个数组。
            FST_meta{n+1} = [FST_meta{n+1}, [repmat(f_meta,1,G_num);G_meta(:,1:G_num)]];
            FST{n+1} = cat(3,FST{n+1},S_next);
        
        end
        
    end
end

