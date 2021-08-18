function [bw,L]=threshold_spheroids_BF(I,varargin)
% this function takes in the BF 10x image and
% identifies spherical objects in the image by adaptive thresholding and
% watershedding. 
% It returns binary and label matrix images.
%optionally it can take in an additional argument which will set the
%threshold value to that argument. 
% author: mnikolic@umd.edu
%% 

t=0.15;
if nargin>1
    t=varargin{1};
end

%%

% I=normalizeImage(I);
bw=imbinarize(I,'adaptive','Sensitivity',0.4);
bw=bwareaopen(bw,200);
bw=imclearborder(bw);
bw=imdilate(bw,strel('disk',3));
bw=imfill(bw,'holes');
bw=imerode(bw,strel('disk',3));
% imagesc(bw);
D = bwdist(~bw);
D = -D;
D = imhmin(D,0.3);
L = watershed(D);
L(~bw) = 0;
% imagesc(L)
% colormap parula
%%
plotON=0;
if plotON
    figure(7)
    subplot(1,2,1);
    imagesc(L);
    axis tight
    axis equal
    subplot(1,2,2);
    B=bwboundaries(L);
    imagesc(I);
    axis tight
    axis equal
    hold on
    for k=1:length(B)
        b=B{k};
        plot(b(:,2),b(:,1),'r-');
    end
    hold off
    drawnow;
end
bw=L>0;
% stats=regionprops(L,scan_image,'MeanIntensity','PixelValues');
%
% mean_shift=[stats.MeanIntensity];
% pixel_shift={stats.PixelValues};
% aggregate_shift=scan_image(L>0);
% Nobjects=length(stats);
end
