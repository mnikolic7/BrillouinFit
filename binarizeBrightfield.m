function [bw] = binarizeBrightfield(I,varargin)
%This function uses stdfilt function from matlab to analyze local
%variations in the brightfield image and then finds the objects. It returns
%the binary image that contains 1s where there are objects in the image. It
%also performs a few morphological operations on the std filt image to
%perform that. dependencies: normalizeImage.m and stdfilt (from matlab).
%   author: mnikolic@umd.edu
%%
I=normalizeImage(double(I));
%%
%
J=stdfilt(I);
%
% imagesc(J)
% colorbar
%
% bw=J>0.02; <- this is used for 20x shape_bf images.
% bw=J>0.045; %<- this is used for cropped 60x images. 
bw=J>0.045;
% T = graythresh( J );
% bw=J>T;
% histogram(J);
% title(num2str(T));
% waitforbuttonpress;
%adjust the above threshold accordingly, 
%but save which threshold was used for what

bw=imdilate(bw,strel('disk',15));
bw=imfill(bw,'holes');
bw=bwareaopen(bw,1500); %previously I had 2000 for the 20x images here.

bw=imerode(bw,strel('disk',15));
bw=imclearborder(bw);
bw = bwareafilt(bw,1);
I2=I;
%%
background=I(~bw);

pd = fitdist(background,'Normal');
f=1.5;
t1=pd.mu-f*pd.sigma;
t2=pd.mu+f*pd.sigma;

bw2=(I>t1).*(I<t2);
bw2=1-bw2; % don't know why this line is here but it makes the thing works so I will keep it.


bw2=bwareaopen(bw2,20);

bw2=imdilate(bw2,strel('disk',8));

bw2=imfill(bw2,'holes');

bw2=bwareaopen(bw2,1500); %previously I had 2000 for the 20x images here.


bw2=imerode(bw2,strel('disk',10));

%%
bw2=imclearborder(bw2);
%%
% imagesc(bw2)
%%
bw2 = bwareafilt(bw2,1);
% imagesc(bw2)
% colorbar
bw=bw2;


%%
% I2(bw2)=1;
if nargin>1
    if strcmpi(varargin{1},'plotOn')
        hf=figure(9);
        hf.Position=[22.3333 209.6667 560 420];
        imagesc([I2,bw])
        colormap gray
        axis equal
    end
end
end

