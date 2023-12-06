function output = eval_layer_1D_FB(input, window, down, freqs, convtype)
% 
% INPUT
% input         array whose columns are assumed to be the input signals
% window        a single window
% convtype      type of convolution
% down          decimation parameter
% freq          frequency sampling set

num_filters = length(freqs);
output_size = length(STFT_1D_FB(input(:,1), window, freqs, convtype, down)); 
output = zeros(output_size, num_filters * size(input,2)); 

% compute the short time Fourier transform for each column of input
for kk = 1:size(input,2)
    outidx = (1:num_filters) + (kk-1)*num_filters;
    output(:,outidx) = STFT_1D_FB(input(:,kk), window, freqs, convtype, down);
end
