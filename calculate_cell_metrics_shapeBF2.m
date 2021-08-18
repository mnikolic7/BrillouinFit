function [metrics_struct] = calculate_cell_metrics_shapeBF2(data_struct)
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
%
% This function also uses the stdfilt on the brightfield image to estimate
% the size and shape of the spheroid
%
%   author: mnikolic@umd.edu


metrics_struct =   struct('name',[],...
                          'mean',[],...
                          'median',[],...
                          'bins',[],...
                          'counts',[],...
                          'FSR',[],...
                          'GHz_per_pixel',[],...
                          'KDE',[],...
                          'min_val',[],...
                          'max_val',[],...
                          'mean_medium',[],...
                          'std_medium',[],...
                          'linewidth',[]);
%%                    
bin_width=0.015; %GHz
bin_centers=6.01:bin_width:6.7;
bin_edges=[bin_centers-bin_width/2 bin_centers(end)+bin_width/2];
%
% W=1;
Good='No';
for W=1:length(data_struct)
    metrics_struct(W).name=data_struct(W).name;
    cellData=data_struct(W).cell;
    metrics_struct(W).mean=nanmean(cellData);
    cellData2=cellData;
    cellData2(cellData2<6.2)=nan;
    metrics_struct(W).mean_T=nanmean(cellData2);
    
    metrics_struct(W).std_shift=nanstd(cellData);
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
    metrics_struct(W).FSR=data_struct(W).FSR;
    metrics_struct(W).GHz_per_pixel=data_struct(W).GHz_per_pixel;
    
    
    params=data_struct(W).params;
    all_widths=data_struct(W).GHz_per_pixel*0.5*(params(2,:)+params(5,:));
    bw_mask=data_struct(W).mask_cell;
    
    scan_image_w=reshape(all_widths, size(data_struct(W).scan_image,1), []); %scan_image is the result. Just plot this to see the result
    scan_image_w=rot90(scan_image_w);
    
    scan_image_w(~bw_mask)=nan;
    metrics_struct(W).linewidth=nanmean(scan_image_w(:));
    %%
%     I=data_struct(W).full_brightfield;
%     rect_crop=data_struct(W).rect_crop;
%     rect_roi=data_struct(W).rect_roi;

    I_cropped=data_struct(W).brightfield;
%     imagesc(I_cropped);
    %%
    bw=binarizeBrightfield(I_cropped,'PlotOn');
    B=bwboundaries(bw);
    hold on
    b=B{1};
    plot(b(:,2),b(:,1),'r-');
    hold off
    Good = questdlg('Is the BF thresholding good?', ...
            'Thresholding good?','Yes','No','No');
%     Good
    %
    S = regionprops(bw,'Area','Centroid','Circularity','MaxFeretProperties');
    A=[S.Area];
    [A,idx]=max(A);
    S=S(idx);
%     Centroid=S(idx).Centroid;
    Circularity=S.Circularity;
    MaxFeretDiameter=S.MaxFeretDiameter;
%     MFC=S(idx).MaxFeretCoordinates;
    %%
%     imagesc(I_cropped);
%     B=bwboundaries(bw);
%     hold on
%     b=B{1};
%     plot(b(:,2),b(:,1),'r-');
%     plot(Centroid(2),Centroid(1),'o');
%     plot(MFC(:,2),MFC(:,1),'g-');
%     hold off
    %%
    if strcmpi(Good,'Yes')
        metrics_struct(W).bf_area=A;
        metrics_struct(W).bf_circ=Circularity;
        metrics_struct(W).bf_MaxFeretDiameter=MaxFeretDiameter;
    else 
        metrics_struct(W).bf_area=nan;
        metrics_struct(W).bf_circ=nan;
        metrics_struct(W).bf_MaxFeretDiameter=nan;
    end
    %% 
end
end

