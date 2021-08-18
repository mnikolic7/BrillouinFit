function [bw,L]=manually_pick_size_cellsd0_BF(I)
% This is a special variation of the function to analyze brightfield images
% taken when imaging cells on day 0.
% this function takes in the BF 10x image and
% identifies spherical objects in the image by adaptive thresholding and
% watershedding. 
% It returns binary and label matrix images.
%optionally it can take in an additional argument which will set the
%threshold value to that argument. 
% author: mnikolic@umd.edu
%% 


%%

L=zeros(size(I));
hf=figure(1);
hf.Position=[21 55.6667 1.2413e+03 574];
imagesc(I);
colormap gray
axis square
for k=1:20
    title(num2str(k));
    h=imfreehand;
    bw_curr=createMask(h);
    L(bw_curr)=k;
end
bw=L>0;
end
