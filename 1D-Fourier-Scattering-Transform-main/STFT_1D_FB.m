function Vf = STFT_1D_FB(input, window, freq, convtype, down)
% computes the short-time Fourier transform of f with window g and using a
% filter bank implementation
%
% INPUTS
% input         samples of input function
% window        window function
% convtype      type of convolution
% down          downsampling parameter
% freq          frequency sampling set

% make in columns
input = input(:);
window = window(:);
freq = freq(:);

M = length(window);
N = length(input);

if nargin < 5
    down = 1;
end

switch convtype
    
    case 'circular'
        % create matrix of shifts
        input = [input; input(1:M-1)];
        
        outputlen = ceil(N/down);
        S = zeros(M,outputlen);
        for kk = 1:outputlen
            startidx = 1+(kk-1)*down;
            S(:,kk) = input(startidx:(startidx+M-1));
        end
        
    case 'nonperiodic'
        % create matrix of shifts
        outputlen = ceil((N-M+1)/down); 
        S = zeros(M,outputlen);
        for kk = 1:outputlen
            startidx = 1+(kk-1)*down;
            S(:,kk) = input(startidx:(startidx+M-1));
        end
        
    otherwise
        error('invalid convolution type')
end

% Wconj = conj(dftmtx(M))/sqrt(M) % not what is expected;
Wdiag = diag(window);
Fmatrix = exp(2*pi*1i*(freq*(0:(M-1)))/M)/sqrt(M);

% after the transpose, this orientation has time on the vertical, filter
% banks on the horizontal
Vf = (Fmatrix * Wdiag * S)';
