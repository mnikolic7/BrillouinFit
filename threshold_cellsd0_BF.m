function [bw,L]=threshold_cellsd0_BF(I)
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

% I=normalizeImage(I);
T=adaptthresh(I);
I2=abs(I-T);
% imagesc(I2)
% axis square
%
bw=edge(I2,'canny',0.1);
% imagesc(bw);
%
bw=imdilate(bw,strel('disk',8));
bw=imfill(bw,'holes');
bw=imerode(bw,strel('disk',8));
bw=imclearborder(bw);
% bw=bwareaopen(bw,200);
bw=imgaussfilt(double(bw),5);
bw=imbinarize(bw);
bw=bwareaopen(bw,300);
% imagesc(bw);
%
D = bwdist(~bw);
D = -D;
D = imhmin(D,0.35);
L = watershed(D);
L(~bw) = 0;
% imagesc(L)
%%
stats=regionprops(L,'Area');
area=[stats.Area];

normArea=area./mean(area);
idx=(normArea<0.6 | normArea>1.5);
area2=area;
area2(idx)=nan;
% plot(area,'.-');
% hold on
% plot(area2,'o')
% hold off
%
L2=L;
for k=1:length(idx)
    if idx(k)
        L2(L2==k)=0;
    end
end
%%
L=L2;
bw=L2>0;
%%
plotON=1;
if plotON
    figure(7)
    subplot(1,2,1);
    imagesc(L);
    axis tight
    axis equal
    ax2=subplot(1,2,2);
    B=bwboundaries(L);
    imagesc(I);
    axis tight
    axis equal
    colormap(ax2,'gray');
    hold on
    for k=1:length(B)
        b=B{k};
        plot(b(:,2),b(:,1),'r-');
    end
    hold off
    drawnow;
end

end
