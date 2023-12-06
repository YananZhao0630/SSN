function [ PC1,PC2 ] = PCALine( input1, input2 )
% Perform the PCA in line (reading the image line by line) 
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Signal Processing Lab, University of Iceland
% 11/10/2014
% ---------------------------------

hdr_in1 = envihdrread([input1,'.hdr']);
hdr_in2 = envihdrread([input2,'.hdr']);


if strcmp(hdr_in1.interleave, 'bip') == 0 || strcmp(hdr_in2.interleave, 'bip') == 0
    disp('The input images have to be in ENVI BIP format!')
    return
end
[precision1, machineformat1] = envInfo(hdr_in1);
[precision2, machineformat2] = envInfo(hdr_in2);


nrow = hdr_in1.lines;
ncol = hdr_in1.samples;
nband = hdr_in1.bands;
size_n = ncol*nband;


[cov_mat,mean_vec] = covmatEval(input1, input2);
cmat1 = cov_mat(1:nband,1:nband);
cmat2 = cov_mat(nband+1:end,nband+1:end);

[coeff1, eigenvalues1] = pcacov(cmat1);
[coeff2, eigenvalues2] = pcacov(cmat2);

avar1 = (mean_vec(1:nband));
avar2 = (mean_vec(nband+1:end));

%new_coeff = coeff(:,1:3);
w_eigen1 = eigenvalues1/sum(eigenvalues1);
w_eigen2 = eigenvalues2/sum(eigenvalues2);

for n1=1: size(coeff1,2)
    %new_coeff(:,i) = coeff(:,i);
    if sum(w_eigen1(1:n1)) > 0.99
        break;
    end
end

for n2=1: size(coeff2,2)
    %new_coeff(:,i) = coeff(:,i);
    if sum(w_eigen2(1:n2)) > 0.99
        break;
    end
end
n = max(n1,n2);

if (n < 5)
    if nband < 5
        n = nband;
    else
        n = 5;
    end
elseif (n > 10) 
    n = 10;
end

new_coeff1 = coeff1(:,1:n);
new_coeff2 = coeff2(:,1:n);

PC1='PC1';
PC2='PC2';

fileIN1 = fopen(input1, 'r');
fileIN2 = fopen(input2, 'r');
filePC1 = fopen(PC1,'w');
filePC2 = fopen(PC2,'w');
for r = 1 : nrow
    %%%%%%%%%%%%%%% the first line of the b-th band
    line1 = fread(fileIN1, size_n, precision1, 0, machineformat1);
    line1 = reshape(line1, nband, ncol);
    
    line2 = fread(fileIN2, size_n, precision2, 0, machineformat2);
    line2 = reshape(line2, nband, ncol);
    
    %%%%%%%%%%%%%%% mean subtraction
    line1 = line1 - repmat(avar1',1,ncol);
    line2 = line2 - repmat(avar2',1,ncol);

    
    %%%%%%%%%%%%%%% principal components
    Pc1 = new_coeff1' * line1;
    Pc2 = new_coeff2' * line2;
    
    % bip
    Pc1 = reshape(Pc1,1,ncol*n);
    fwrite(filePC1,Pc1,class(Pc1));
    
    Pc2 = reshape(Pc2,1,ncol*n);
    fwrite(filePC2,Pc2,class(Pc2));
end
fclose(fileIN1);
fclose(fileIN2);
fclose(filePC1);
fclose(filePC2);

hdrWrite(PC1,nrow,ncol,n,class(Pc1));
hdrWrite(PC2,nrow,ncol,n,class(Pc2));

end

