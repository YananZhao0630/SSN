function G = mygausswin(len)
% generates a vector containing uniform samples of a standard Gaussian on
% the interval [-2,2] normalized to have max of 1
% INPUT 
% len   length of the desired window

t = linspace(-2,2,len); 
G = exp(-t.^2/2); 
G = G(:)/norm(G);
