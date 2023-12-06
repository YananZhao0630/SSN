function Z = ChiSquareLine(input,filename)

% function for the computation of the Chi square distribution based on EM algorithm for Gaussian mixtures
% estimation
%
%----------------------------------
% input:
% - input       - string of the whole path of the ENVI bip format image
% - filename
% 
% output: 
% - Z           - Chi square distribution of the data 
%
%----------------------------------
% Nicola Falco 
% nicolafalco@ieee.org
% 
% Prashanth Reddy Marpu
% prashanthmarpu@ieee.org
% 
% Signal Processing Lab, University of Iceland
% 20/11/2011 first version
% 11/10/2015 last revision

% ---------------------------------

hdr_in = envihdrread([input,'.hdr']);

if strcmp(hdr_in.interleave, 'bip') == 0
    disp('The input images have to be in ENVI BIP format!')
    return
end

[precision, machineformat] = envInfo(hdr_in);

nrow = hdr_in.lines;
ncol = hdr_in.samples;
nband = hdr_in.bands;
size_n = ncol*nband;

dim = min(50000,nrow*ncol);
v(:,1) = randi([1 nrow],1,dim);
v(:,2) = randi([1 ncol],1,dim);
[idx,IX] = sort(v(:,1));
new = [v(IX,1),v(IX,2)];

for b = 1 : nband
    fprintf('\nBand %d: ',b);
    file_in = fopen(input,'r');
    
    rand_data = [];
    for r = 1 : nrow
        line = fread(file_in, size_n, precision, 0, machineformat);
        line = reshape(line, nband, ncol);
        
        rand_data = [rand_data ,line(b,new(idx==r,2))];
    end
    fclose(file_in);
    
    % EM_GM Algorithm %
    disp('EM_GM Algorithm in progress . . . ');
    [W,M,V] = EM_GM(rand_data(:),3,[],2,[],[]);
    [n,indi] = min(abs(M));
    std(b) = sqrt(V(indi));
end
fprintf('\n');

file_in = fopen(input,'r');
file_out = fopen(filename,'w');
for r = 1 : nrow
    line = fread(file_in, size_n, precision, 0, machineformat);
    line = reshape(line, nband, ncol);
    
    Z = (line./ repmat(std',1,ncol)) .^2;
    
    Z = sum(Z,1);
    fwrite(file_out,Z,class(Z));
end
fclose(file_in);
fclose(file_out);

hdrWrite(filename,nrow,ncol,1,class(Z));
clear Z
Z = enviread(filename);
save(filename, 'Z');
end