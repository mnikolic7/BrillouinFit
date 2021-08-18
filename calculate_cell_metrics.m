function [metrics_struct] = calculate_cell_metrics(data_struct,varargin)
%CALCULATE_CELL_METRICS takes the data_struct from a condition and
%measurement and extracts cell data from the selected cell, and calculates
%parameters that can be used to compare with other conditions.
%parameters/metrics are listed as fields of the output struct.
%   author: mnikolic@umd.edu

%% handle variable arguments in: analysis options
select_manual_BF=0;
select_auto_BF=0;
select_auto_Brillouin=0;
good='No';
if nargin==2
    select_manual_BF=varargin{1};
elseif nargin==3
    select_manual_BF=varargin{1};
    select_auto_BF=varargin{2};
elseif nargin==4
    select_manual_BF=varargin{1};
    select_auto_BF=varargin{2};
    select_auto_Brillouin=varargin{3};
else
    %     disp 'extra arguments (after 4th) will be ignored';
end

%%
%initialize output struct
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
%% initialize variables for fitting histogram
bin_width=0.015; %GHz
bin_centers=6:bin_width:6.8;
bin_edges=[bin_centers-bin_width/2 bin_centers(end)+bin_width/2];
%% for each entry in the data_struct calculate the metrics and save
for W=1:length(data_struct)
    %save metadata like name and data path
    metrics_struct(W).name=data_struct(W).name;
    
    %% list of shifts, mean, and std in the cell.
    %do three different ways:
    %   1. just the handpicked average of all shifts from main (it removes
    %   all shifts larger than mean(medium)+3.5*std(medium) and does no
    %   image manipulation. So whatever is >threshold is included.
    %   2. global threshold 6.15 GHz
    %   3. global threshold 6.4 GHz
    %extract a list of shifts inside the cell
    cellData=data_struct(W).cell;
    metrics_struct(W).cellData=cellData;
    %extract a list of all shifts inside the brillouin scan image
    cellDataFull=data_struct(W).scan_image;
    cellDataFull=cellDataFull(:);
    metrics_struct(W).cellDataFull=cellDataFull;
    %save the average of the automatically thresholded cell image.
    %mean of all shifts that are 3.5sigma>medium shifts.
    metrics_struct(W).mean=nanmean(cellData);
    
    scan_image=data_struct(W).scan_image;
    x_scan=data_struct(W).x_scan;
    y_scan=data_struct(W).y_scan;
    
    
    metrics_struct(W).x_scan=x_scan;
    metrics_struct(W).y_scan=y_scan;
    metrics_struct(W).scan_image=scan_image;
    
    
    bw6p4=scan_image>6.4;
    bw6p4=imdilate(bw6p4,strel('disk',1));
    bw6p4=imfill(bw6p4,'holes');
    bw6p4=imerode(bw6p4,strel('disk',1));
    
    bw6p15=scan_image>6.15;
    bw6p15=imdilate(bw6p15,strel('disk',1));
    bw6p15=imfill(bw6p15,'holes');
    bw6p15=imerode(bw6p15,strel('disk',1));
    
    bw6p3=scan_image>6.3;
    
    %save all shifts with global thresholds
    cellData2=scan_image(bw6p15);
    cellData3=scan_image(bw6p4);
    cellData4=scan_image(bw6p3);
    
    %save means of fixed threshold shifts.
    metrics_struct(W).mean_T6p15=nanmean(cellData2(:));
    metrics_struct(W).mean_T6p4=nanmean(cellData3(:));
    metrics_struct(W).mean_T6p3=nanmean(cellData4(:));  
    %note 6p3 threshold is meant to ignore holes and will not worry about
    %where the >6p3 pixels are or whether there is a hole in the image
    %(lumen). The 6p15 and 6p4 shifts use hole filling, so if there is a
    %hole in the image the average will include that hole as well. 
    
    %save standard deviation of all shifts in the cell, also the fixed
    %threshold ones
    metrics_struct(W).std_shift=nanstd(cellData);
    metrics_struct(W).std_shift_T6p15=nanstd(cellData2(:));
    metrics_struct(W).std_shift_T6p4=nanstd(cellData3(:));
    
    %save the median
    metrics_struct(W).median=nanmedian(cellData);
    %save the min and max value
    metrics_struct(W).min_val=min(cellData);
    metrics_struct(W).max_val=max(cellData);
    %% medium metrics from initial selection
    medium=data_struct(W).medium;
    metrics_struct(W).mean_medium=nanmean(medium);
    metrics_struct(W).std_medium=nanstd(medium);
    metrics_struct(W).FSR=data_struct(W).FSR;
    metrics_struct(W).GHz_per_pixel=data_struct(W).GHz_per_pixel;
    
    %% histogram and KDE
    %save the bins
    metrics_struct(W).bins=bin_centers;
    %save histogram
    [counts,~]=histcounts(cellData,bin_edges,'Normalization','probability');
    metrics_struct(W).counts=counts;
    %calculate and save KDE
    pdKDE = fitdist(cellData,'kernel','BandWidth',bin_width);
    yKDE=pdf(pdKDE, bin_centers)*bin_width;
    metrics_struct(W).KDE=yKDE;
    
    %save histogram and KDE for full image as well (including medium)
    [countsFull,~]=histcounts(cellDataFull,bin_edges,'Normalization','probability');
    metrics_struct(W).countsFull=countsFull;
    pdKDEFull = fitdist(cellDataFull,'kernel','BandWidth',bin_width);
    yKDEFull=pdf(pdKDEFull, bin_centers)*bin_width;
    metrics_struct(W).KDEFull=yKDEFull;
    
    %note: if you take the fit of the medium distribution,normalize it properly,
    % and subtract it from the full histogram/KDE you can get another
    % estimate for the shift distribution.
    
    %% save all parameters from the fitting.
    % include mean linewidth.
    params=data_struct(W).params;
    all_widths=data_struct(W).GHz_per_pixel*0.5*(params(2,:)+params(5,:));
    bw_mask=data_struct(W).mask_cell;
    
    scan_image_w=reshape(all_widths, size(data_struct(W).scan_image,1), []); %scan_image is the result. Just plot this to see the result
    scan_image_w=rot90(scan_image_w);
    
    scan_image_w0=scan_image_w;
    scan_image_w2=scan_image_w;
    scan_image_w3=scan_image_w;
    scan_image_w4=scan_image_w;
    scan_image_w2(~bw6p15)=nan;
    scan_image_w3(~bw6p4)=nan;
    scan_image_w0(~bw_mask)=nan;
    scan_image_w4(~bw6p3)=nan;
    metrics_struct(W).linewidth=nanmean(scan_image_w0(:));
    metrics_struct(W).linewidth6p15=nanmean(scan_image_w2(:));
    metrics_struct(W).linewidth6p4=nanmean(scan_image_w3(:));
    metrics_struct(W).linewidth6p3=nanmean(scan_image_w4(:));

    
    %% plot_radial_profile_of_brillouin shifts
    %first ask user if this image is a good one. 
    W
    rect_roi=data_struct(W).rect_roi;
    if (rect_roi(3)==1 || rect_roi(4)==1)
        disp(['omitting ', data_struct(W).name,' because it is xz or yz']);
        metrics_struct(W).radial_profile=nan(200,1);
        metrics_struct(W).tan_profile=nan(200,1);
    else
        %
        [radial_image,r_scan,th_scan] = plot_radial_profile(x_scan,y_scan,scan_image,bw6p15);
        radial_profile=nanmean(radial_image,2);
        radial_image2=radial_image;
        radial_image2(radial_image2<6.4)=nan;
        tan_profile=nanmean(radial_image2,1);
        
        metrics_struct(W).radial_profile=radial_profile;
        metrics_struct(W).tan_profile=tan_profile;
        metrics_struct(W).radial_image=radial_image;
        %     metrics_struct(W).radial_profile_err=nanstd(radial_image,2);
        %     metrics_struct(W).tan_profile_err=nanstd(radial_image,1);
        metrics_struct(W).r_scan=r_scan;
        metrics_struct(W).th_scan=th_scan;
    end
    %% (optional) manually select shape in BF image
    if select_manual_BF==1
        
        name=data_struct(W).name;
        
        full_brightfield=data_struct(W).full_brightfield;
        rect_crop=data_struct(W).rect_crop;
        rect_roi=data_struct(W).rect_roi;
        %
        
        %plot brillouin so that you can double verify it is the right cell
        hf6=figure(6);
        hf6.Position=[17.6667  383.6667  267.3333  252.6667];
        ax=axes(hf6);
        imagesc(ax,scan_image);
        colormap jet
        caxis([6 6.5]);
        colorbar;
        axis equal
        bw_b=data_struct(W).mask_cell;
        bw_b=bwareafilt(bw_b,1);
        BB=bwboundaries(bw_b);
        bb=BB{1};
        hold on
        plot(ax,bb(:,2),bb(:,1),'r-','linewidth',1);
        hold off
        
        
        %% ask user to select Brightfield image for shape identification.
        rect_roi_full=rect_roi;
        rect_roi_full(1)=rect_roi_full(1)+rect_crop(1);
        rect_roi_full(2)=rect_roi_full(2)+rect_crop(2);
        
        axlim=[rect_crop(2)-100,rect_crop(2)+rect_crop(4)+200, rect_crop(1)-100,rect_crop(1)+rect_crop(3)+200];
        
        mask_bf=detect_cells_in_brightfield(full_brightfield,axlim,rect_roi_full,name);
        %% calculate shape data
        pixel_size=1; %change based on the objective.
        stats = regionprops(mask_bf,'Area','MajorAxisLength','MinorAxisLength','Perimeter','BoundingBox');
        
        areas=[stats.Area];
        [~,idx]=max(areas);
        
        area=stats(idx).Area*pixel_size.^2;
        aspect_ratio=stats(idx).MajorAxisLength./stats(idx).MinorAxisLength;
        circularity=(4*area*pi)/(stats(idx).Perimeter*pixel_size).^2;
        major_axis=stats(idx).MajorAxisLength;
        mask_bf_cropped=imcrop(mask_bf,stats(idx).BoundingBox);
        
        
        %brightfield shape measurements
        metrics_struct(W).area=area;
        metrics_struct(W).aspect_ratio=aspect_ratio;
        metrics_struct(W).circularity=circularity;
        metrics_struct(W).major_axis=major_axis;
        metrics_struct(W).mask_bf=mask_bf_cropped;
        close(hf6);
    end
    %% (optional) automatically select shape in BF image
    if select_auto_BF==1
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
        %         imagesc(I_cropped);
        %         B=bwboundaries(bw);
        %         hold on
        %         b=B{1};
        %         plot(b(:,2),b(:,1),'r-');
        %         plot(Centroid(2),Centroid(1),'o');
        %         plot(MFC(:,2),MFC(:,1),'g-');
        %         hold off
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
    end
    %% (optional) automatically select shape in Brillouin image
    if select_auto_Brillouin==1
        x_scan=data_struct(W).x_scan;
        y_scan=data_struct(W).y_scan;
        
        xumpp=x_scan(2)-x_scan(1);
        yumpp=y_scan(2)-y_scan(1);
        bw_mask=data_struct(W).mask_cell;
        bw=bwareaopen(bw_mask,25);
        CC=bwconncomp(bw);
        S=regionprops(CC,'Area','Perimeter','MaxFeretProperties','MinFeretProperties');
        % keep only largest area in image in case there are outliers;
 
        
        [Area,iiidx]=max([S.Area]);
        Perimeter=max([S.Perimeter]);
        MaxFeretDiameter=S(iiidx).MaxFeretDiameter;
        MinFeretDiameter=S(iiidx).MinFeretDiameter;
        
        
        metrics_struct(W).BArea=Area*(xumpp)^2;
        metrics_struct(W).BPerimeter=Perimeter*xumpp;
        metrics_struct(W).BCirc=(Perimeter.^2)/(4*pi*Area);
        metrics_struct(W).BDiameter=0.5*(MaxFeretDiameter+MinFeretDiameter);
    end
end
end

