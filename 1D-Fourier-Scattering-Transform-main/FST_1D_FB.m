function [S, U, Smeta] = FST_1D_FB(input, depth, windows, dec, freqs, convtype)
% calculates 1D FST using a filter bank implementation
% input     array of samples of function
% depth     number of layers in scattering
% win       a cell of filters for each layer and the final filtering
%           should be at least depth + 1
%           the d+1 element is the filter used for the final convolution
% dec       array of length d+1 specifies how much to downsample the final
%           coefficients
% freqs     cell corresponding to which frequencies to sample
%           controls the number of filters used per layer
% convtype  type of convolution used in each layer


input = input(:);

% coefficients for the zeroth order scattering
U{1} = input; 

% zeroth order scattering is just a smoothing
S{1} = STFT_1D_FB(input, windows{1}, 0, convtype,dec(1)); 
Smeta{1} = [];

% iterate over each layer
for d = 1:depth
    
    % window for the d-th layer
    g = windows{d};
    
    % number of filters in the d-th layer
    numbanks = length(freqs{d});
    
    % number of coefficients in the previous layer
    curr_num_signals = size(U{d},2);
        
    % update the meta cell
    if d == 1
        Smeta{d+1} = reshape(freqs{d},1,[]);    % row vector of frequencies 
    else 
        Smeta{d+1} = reshape(repmat(Smeta{d},numbanks,1),d-1,[]); 
        size(Smeta{d+1});
        Smeta{d+1} = [Smeta{d+1}; repmat(freqs{d},1,curr_num_signals)];
    end
    
    % intermediate coefficients, no spatial downsampling
    U{d+1} = abs(eval_layer_1D_FB(U{d}, g, 1, freqs{d}, convtype));
    
    % scattering coefficients possibly with spatial downsampling 
    S{d+1} = eval_layer_1D_FB(U{d+1}, g, dec(d), 0, convtype); 
end

% now do all the low pass
%[Phi, S, Smetadata] = U_to_features(U, Umetadata, depth, win{depth+1}, dec(depth+1));
