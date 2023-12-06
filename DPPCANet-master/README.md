# DPPCANet



## Introduction

DPPCANet is a robust deep learning method for change detection of unbalanced multi-temporal synthetic aperture radar images, which mainly includes 1) generating difference maps; 2) parallel FCM clustering to provide training sample pseudo-labels; 3) constructing over-sampling and under-sampling Pixel classification based on sampled PCANet + SVM model.



## Requirements

MATLAB 2018a 



## function

Weighted-pooling convolution:

`I_wp = WP(I,k)`  



Generate log-ratio image:

`I_lr = di_gen(I_1,I_2)`



Cumulative weighted pooling:  

T is the Cumulative times

`M = Normalized(matrix)` is a robust normalization method. The average value of the 50 elements with the highest value in the input matrix is the upper bound, and the average value of the 50 elements with the lowest value is the lower bound. The matrix is soft normalized.

`DDI = Normalized(CWP(I_ori,T))`



Gabor wavelet transform feature extraction:

`[f1,f_1] = Gabor_fea(I_map)`



Parallel clustering: two sets of mapped DDI I_map1,I_map2 and Gabor feature vectors f_1 and f_2

`im_lab = paralleclustering(f_1,I_map1,f_2,I_map2)`



FCM:

`[center, U, obj_fcn] = FCM(data, cluster_n, options)`

<!--data ---- n * m matrix, representing n samples, each sample has m-dimensional eigenvalues-->
<!--cluster_n ---- scalar, representing the number of aggregation centers, that is, the number of categories-->
<!--options ---- 4 * 1 column vector, where-->
<!--options (1): exponent of membership matrix U,> 1 (default: 2.0)-->
<!--options (2): maximum number of iterations (default: 100)-->
<!--options (3): minimum change in membership, iteration termination condition (default: 1e-5)-->
<!--options (4): whether to output information flags for each iteration (default: 0)-->

<!--center ---- cluster center-->
<!--U ---- membership matrix-->
<!--obj_fcn ---- objective function value-->



The train and test of PCANet+SVM:

`PCANet_SVM_train;`

`PCANet_SVM_test;`



Performance evaluation:

`[TP,TN,FP,FN,MC,MU,FPR,FNR,OER,PCC,Kappa] = PE(im_res,im_gt)`



## note:

Because the training samples are randomly selected, the results obtained by each run are not fixed, and the results obtained are slightly different.



## Acknowledgements:

Thanks to F.Gao for the basic code about PCANet, and it has inspiration for DPPCANet.



















