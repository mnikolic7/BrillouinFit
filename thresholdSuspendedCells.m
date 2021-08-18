function [mean_shift,pixel_shift,aggregate_shift,Nobjects]=thresholdSuspendedCells(scan_image,varargin)
% this function takes in the brillouin shift image (scan_image) and
% identifies spherical objects in the image by thresholding and
% watershedding. 
% It returns mean_shift - a vector of average shift in each image
%            pixel_shift - a cell with pixel shift values for each
%            object/cell
%            aggregate_shift - a vector of all shifts in the image after
%            thresholding
%optionally it can take in an additional argument which will set the
%threshold value to that argument. 
% author: mnikolic@umd.edu
%% 

t=6.15;
if nargin>1
    t=varargin{1};
end

%%

I=scan_image;
% I=normalizeImage(I);
bw=imbinarize(I,t);
bw=bwareaopen(bw,5);
% imagesc(bw);
D = bwdist(~bw);
D = -D;
L = watershed(D);
L(~bw) = 0;
% imagesc(L)
%
plotON=0;
if plotON
    figure(7)
    subplot(1,2,1);
    imagesc(L);
    axis tight
    axis equal
    subplot(1,2,2);
    B=bwboundaries(L);
    imagesc(scan_image);
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
stats=regionprops(L,scan_image,'MeanIntensity','PixelValues');
%
mean_shift=[stats.MeanIntensity];
pixel_shift={stats.PixelValues};
aggregate_shift=scan_image(L>0);
Nobjects=length(stats);
end
