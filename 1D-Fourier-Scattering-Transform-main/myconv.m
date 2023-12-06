function F = myconv(f,g,par)
% computes the convolution 
% INPUT
% f         function
% g         window
% par       type of convolution
%
% Weilin, October 2017

switch par
    case 'circular'
        F = cconv(f,g,length(f));
    case 'same'
        F = conv(f,g,par);
    case 'valid'
        F = conv(f,g,par);
    otherwise
        error('invalid convolution type')
end