% input
addpath('./Scattering_v2');
addpath('./Utilities');
clear all; close all; clc;
s1 = 2^2;  u1 = 10;  
s2 = 2^2;  u2 = 10;  
v1 = 4*pi/s1; v2 = 4*pi/s2;
t1 = 0:1/(2*v1):5;  t1_1 = (1)*t1;
t2 = 0:1/(2*v2):5;  t2_2 = (1)*t2;
T1 = ((t1_1-u1)./s1);
T2 = ((t2_2-u2)./s2);
g1 = 2.^(8).*exp(- pi.*(T1.^(2)));
g_r1 = (1./(sqrt(s1))).*g1.*exp(1j.*v1.*t1_1);
g2 = 2.^(8).*exp(- pi.*(T2.^(2)));
g_r2 = (1./(sqrt(s2))).*g2.*exp(1j.*v2.*t2_2);
G1 = real(g_r1'*g_r2);
N = size(G1,1);
radius = pi/8;
scales = 8;
depth = 1;
%[G0,Gp,G_meta,G_radius] = Gabor_Bank_2D(N,radius,scales);
dec = 'yes';
[FST,FST_meta] = Fourier_Scattering_2D(G1,radius,scales,depth,dec);







