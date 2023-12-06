function [S_thres,M_thres] = threshold_coeff(Scat,Meta,cutoff)
% THRESHOLD_COEFF discards coefficients whose norms are smaller than cutoff
%
% INPUT:
%  Scat             (cell) of scattering coefficients
%  Meta             (cell) of information on Scat
%  cutoff           (real) thresholding parameter
%
% OUTPUT:
%  S_thres		    (cell) of thresholded scattering coefficients
%  M_thres          (cell) of information on S_thres
%--------------------------------------------------------------------------
% Weilin Li ~ June 2017

% initialize
S_thres = [];
M_thres = cell(1,length(Scat));
index = 1;

% keep large coefficients
for j = 1:length(Scat)
    for k = 1:size(Scat{j},3)
        energy = norm(Scat{j}(:,:,k),'fro');
        if energy > cutoff
            S_thres(:,:,index) = Scat{j}(:,:,k);
            if j ==1
                M_thres{j} = [];
            else
                M_thres{j} = [M_thres{j}, Meta{j}(:,k)];
            end
            index = index+1;
        end
    end
end

    

    