function [im_nor] = Normalized(im_or)
% Normalize the average of the 50 highest pixel values and 50 lowest pixel values to the upper and lower limits
[x,y] = size(im_or);
im_or_v = reshape(im_or,1,x*y);
[im,index] = sort(im_or_v,'DESCEND');
im_max = sum(sum(im(1,1:100)))/100;
im_min = sum(sum(im(1,x*y-100+1:x*y)))/100;
im_nor = (im_or-im_min)/(im_max-im_min);
end

