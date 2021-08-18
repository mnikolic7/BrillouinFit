function [metrics_struct] = calculate_cell_metrics_shapeBF(data_struct)
%CALCULATE_CELL_METRICS takes the data_struct from a condition and
%measurement and extracts cell data from the selected cell, and calculates
%a few simple parameters that can be used to compare with other conditions.
%parameters are: 
%                           'mean',[],...
%                           'median',[],...
%                           'bins',[],...
%                           'counts',[],...
%                           'KDE',[],...
%                           'min_val',[],...
%                           'max_val',[]);
%   author: mnikolic@umd.edu


metrics_struct =   struct('mean',[],...
                          'median',[],...
                          'bins',[],...
                          'counts',[],...
                          'KDE',[],...
                          'min_val',[],...
                          'max_val',[],...
                          'mean_medium',[],...
                          'area',[],...
                          'aspect_ratio',[],...
                          'circularity',[],...
                          'std_medium',[]);
                    
bin_width=0.015; %GHz
bin_centers=6.01:bin_width:6.7;
bin_edges=[bin_centers-bin_width/2 bin_centers(end)+bin_width/2];
for W=1:length(data_struct)
    cellData=data_struct(W).cell;
    metrics_struct(W).mean=nanmean(cellData);
    metrics_struct(W).median=nanmedian(cellData);
    metrics_struct(W).bins=bin_centers;
    
    [counts,~]=histcounts(cellData,bin_edges,'Normalization','probability');
    metrics_struct(W).counts=counts;
    pdKDE = fitdist(cellData,'kernel','BandWidth',bin_width);
    yKDE=pdf(pdKDE, bin_centers)*bin_width;
    metrics_struct(W).KDE=yKDE;
    metrics_struct(W).min_val=min(cellData);
    metrics_struct(W).max_val=max(cellData);
    
    medium=data_struct(W).medium;
    metrics_struct(W).mean_medium=nanmean(medium);
    metrics_struct(W).std_medium=nanstd(medium);
    %%
    Ibf=data_struct(W).full_brightfield;
    crop_roi=data_struct(W).rect_crop;
    rect_roi=data_struct(W).rect_roi;
    hf=figure(8);
    ax=axes();
    imagesc(Ibf);
    colormap(ax,'gray');
    set(ax,'ylim',[crop_roi(1)-100,crop_roi(1)+crop_roi(3)+200]);
    set(ax,'xlim',[crop_roi(2)-100,crop_roi(2)+crop_roi(4)+200]);
    abs_roi=rect_roi;
    abs_roi(1)=abs_roi(1)+crop_roi(1);
    abs_roi(2)=abs_roi(2)+crop_roi(2);
    rectangle('Position',abs_roi);
    title('please select the cell');
    %
%     crop_rect=getrect();
%     I2=imcrop(Ibf,crop_roi);
    %%
    h=imfreehand;
    bw_bf=createMask(h);
    bw_bf=bwareafilt(bw_bf,1);
    bw_bf=imgaussfilt(double(bw_bf),5);
    bw_bf=imbinarize(bw_bf);
%     %%
% %     I2=Ibf;
% %     I2(~bw_bf)=0;
% %     bw_outline=bwperim(bw_bf);
%     [~,threshold] = edge(I2,'canny');
%     fudgeFactor = 2.5;
%     BWs = edge(I2,'canny',threshold * fudgeFactor);
% %     BWs=BWs-bw_outline;
% %     BWs=imgaussfilt(double(BWs),2);
% %     bw_out=imgaussfilt(double(bw_outline),2);
% %     bw_out=imbinarize(bw_out);
% %     BWs=imbinarize(BWs);
% 
%     imagesc(BWs);
%     colormap('gray');
% %     axis([1440 1620 1350 1650])
%     %
%     se90 = strel('line',5,90);
%     se0 = strel('line',5,0);
%     BWsdil = imdilate(BWs,[se90 se0]);
% 
%     %
%     BWdfill = imfill(BWsdil,'holes');
%     %
%     seD = strel('diamond',2);
%     BWfinal = imerode(BWdfill,seD);
%     BWfinal = imerode(BWfinal,seD);
%     BWfinal=bwareafilt(BWfinal,1);
%     BWfinal=imgaussfilt(double(BWfinal),2);
%     BWfinal=imbinarize(BWfinal);
%     bw_bf=BWfinal;
%     Idisp=I2;
%     bwoutline=bwperim(bw_bf);
%     Idisp(bwoutline)=1;
%     
%     imagesc(Idisp);
%     title('select regions that are to be removed');
%     h=imfreehand;
%     bw_r=createMask(h);
%     waitforbuttonpress;
%     bw_bf(bw_r)=0;
%     bwoutline=bwperim(bw_bf);
%     Idisp=I2;
%     Idisp(bwoutline)=1;
%     imagesc(Idisp);
%     title('click to continue');
%     bw_bf=bwareafilt(bw_bf,1);
%     waitforbuttonpress;
    pixel_size=1; %change based on the objective. 
    stats = regionprops(bw_bf,'Area','MajorAxisLength','MinorAxisLength','Perimeter');
    area=stats.Area*pixel_size.^2;
    aspect_ratio=stats.MajorAxisLength./stats.MinorAxisLength;
    circularity=(4*area*pi)/(stats.Perimeter*pixel_size).^2;
    
    metrics_struct(W).area=area;
    metrics_struct(W).aspect_ratio=aspect_ratio;
    metrics_struct(W).circularity=circularity;
    metrics_struct(W).mask_brightfield=bw_bf;
%     metrics_struct(W).crop_rect=crop_rect;
    metrics_struct(W).brightfield=Ibf;
%%
end
close(hf);
end

