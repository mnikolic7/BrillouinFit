function [metrics_struct] = calculate_cell_metrics_Bshape(data_struct)
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
bin_centers=6:bin_width:6.7;
bin_edges=[bin_centers-bin_width/2 bin_centers(end)+bin_width/2];
Nbins=length(bin_centers);
for W=1:length(data_struct)
    metrics_struct(W).name=data_struct(W).name;
    cellData=data_struct(W).cell;
    metrics_struct(W).mean=nanmean(cellData);
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
    
    
    x_scan=data_struct(W).x_scan;
    y_scan=data_struct(W).y_scan;
    
    xumpp=x_scan(2)-x_scan(1);
    yumpp=y_scan(2)-y_scan(1);
    
    bw=bwareaopen(bw_mask,25);
    CC=bwconncomp(bw);
    S=regionprops(CC,'Area','Perimeter');
    
    Area=max([S.Area]);
    Perimeter=max([S.Perimeter]);
    
    metrics_struct(W).BArea=Area*(xumpp)^2;
    metrics_struct(W).BPerimeter=Perimeter*xumpp;
    metrics_struct(W).BCirc=(Perimeter.^2)/(4*pi*Area);
    
    
        
  
    
    
    %%
%     hf=figure(4);
%     color_indices=round(linspace(1,256,Nbins));
%     cmap=colormap('jet');
%     ax1=subplot(2,1,1);
%     imagesc(data_struct(W).scan_image)
%     colormap jet
%     caxis([6.0 6.35]);
%     axis square
%     ax2=subplot(2,1,2);
%     hb=bar(bin_centers,yKDE,'FaceColor','flat');
%     for k = 1:Nbins
%         hb.CData(k,:) = cmap(color_indices(k),:);
%     end
%     ylim(ax2,[0 1.2*max(counts)]);
%     hold(ax2,'off');
%     waitforbuttonpress;
    %%
%     bw=data_struct.mask_cell;
%     %%
%     bw=bwareaopen(bw,20);
%     B=bwboundaries(bw);
%     imagesc(bw);
%     hold on
%     b=B{1};
%     plot(b(:,2),b(:,1),'r-');
%     hold off
    

end


end

