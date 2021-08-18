function [mask_cell,mask_medium] = calculate_display(x_scan,y_scan,scan_image,rect_roi,mask_cell,mask_medium,name)
%This function creates figures and plots the Brillouin and Brightfield images
%from the data provided by calculate_cell_metrics_v2.m
%it returns the hf a cell array of two figure handles (for brilloui and
% brightfield images, and also ax - cell of two axes handles)
%   author: mnikolic@umd.edu


%% plot Brillouin

%
mask_cell=BrillouinThreshold(scan_image,mask_cell,mask_medium);

scan_type='xy';
if rect_roi(3)<1
    scan_type='yz';
elseif rect_roi(4)<2
    scan_type='xz';
end

B=bwboundaries(mask_cell);
if length(B)>1
    len=zeros(length(B),1);
    for k=1:length(B)
        len(k)=size(B{k},1);
    end
    [~,idx_b]=max(len);
    B=B(idx_b);
end
b=B{1};

Bm=bwboundaries(mask_medium);
if length(Bm)>1
    len=zeros(length(Bm),1);
    for k=1:length(Bm)
        len(k)=size(Bm{k},1);
    end
    [~,idx_b]=max(len);
    Bm=Bm(idx_b);
end
bm=Bm{1};


hf1=figure(1);
ax1=axes;
hf1.Position=[32 230 500 400]; %[x0, y0, width, height]
title(['name: ', name, ' Brillouin']);
h_Brill_img=imagesc(x_scan,y_scan,scan_image);
set(ax1,'ydir','normal');
if strcmpi(scan_type,'xz') || strcmpi(scan_type,'yz')
    set(ax1,'ydir','reverse');
end
hold on
plot(x_scan(b(:,2)),y_scan(b(:,1)),'-','color',[1	0.462745098039216	0.290196078431373],'linewidth',1);
plot(x_scan(bm(:,2)),y_scan(bm(:,1)),'--','color',[1	0.462745098039216	0.290196078431373],'linewidth',1);
axis equal
axis tight
shading flat
colormap(ax1, 'jet');
colorbar();
caxis([6 6.7]);
xlabel([num2str(max(x_scan),'%.0d'), ' \mum']);
ylabel([num2str(max(y_scan),'%.0d'), ' \mum']);
set(ax1,'FontSize',11);
pbaspect([max(x_scan),max(y_scan), 1])
hold off
%% residual code from picking the figure position
% h1=figure(1);
% h1.Position=[32 230 500 400];
% h2=figure(2);
% h2.Position=[550 230 500 400];
%
% h3=figure(3);
% h3.Position=[32 42 500 150];
% set(h3, 'MenuBar', 'none');
% set(h3, 'ToolBar', 'none');

%% plot brillouin histogram
hf3=figure(3);
hf3.Position=[32 42 500 150];
set(hf3, 'MenuBar', 'none');
set(hf3, 'ToolBar', 'none');
ax3=axes(hf3);
bin_width=0.015; %GHz
bin_centers=6.01:bin_width:6.7;
bin_edges=[bin_centers-bin_width/2 bin_centers(end)+bin_width/2];
[counts,~]=histcounts(scan_image(mask_cell),bin_edges);%,'Normalization','probability');
[countsM,~]=histcounts(scan_image(mask_medium),bin_edges);%,'Normalization','probability');

plot(ax3,bin_centers,counts,'-');
hold(ax3,'on');
plot(ax3,bin_centers,countsM,'--');
hold(ax3,'off');

%% select Brillouin image

answer='Yes';
while ~strcmpi(answer,'No')
    answer = questdlg('Do you want to select again?', ...
        'Brillouin Selection Menu', ...
        'Cell manual','Both','No','No');
    % Handle response
    switch answer
        case 'Cell manual'
            h=imfreehand(ax1);
            mask_cell=createMask(h,h_Brill_img);
%             mask_cell=BrillouinThreshold(scan_image,mask_cell,mask_medium);
        case 'Both'
            disp('select cell')
            h=imfreehand(ax1);
            mask_cell=createMask(h,h_Brill_img);
            
            disp('select medium')
            h=imfreehand(ax1);
            mask_medium=createMask(h,h_Brill_img);
            %threshold out by 3.5 sigma
            mask_cell=BrillouinThreshold(scan_image,mask_cell,mask_medium);
        case 'No'
    end
    
    %     cla(ax1);
    B=bwboundaries(mask_cell);
    if length(B)>1
        len=zeros(length(B),1);
        for k=1:length(B)
            len(k)=size(B{k},1);
        end
        [~,idx_b]=max(len);
        B=B(idx_b);
    end
    b=B{1};
    
    Bm=bwboundaries(mask_medium);
    if length(Bm)>1
        len=zeros(length(Bm),1);
        for k=1:length(Bm)
            len(k)=size(Bm{k},1);
        end
        [~,idx_b]=max(len);
        Bm=Bm(idx_b);
    end
    bm=Bm{1};
    
    h_Brill_img=imagesc(ax1,x_scan,y_scan,scan_image);
    hold(ax1,'on');
    plot(ax1,x_scan(b(:,2)),y_scan(b(:,1)),'-','color',[1   0.462745098039216   0.290196078431373],'linewidth',1);
    plot(ax1,x_scan(bm(:,2)),y_scan(bm(:,1)),'--','color',[1    0.462745098039216   0.290196078431373],'linewidth',1);
    %     axis equal
    %     axis tight
    %     shading flat
    colormap(ax1, 'jet');
    colorbar(ax1);
    caxis(ax1,[6 6.7]);
    %     xlabel([num2str(max(x_scan),'%.0d'), ' \mum']);
    %     ylabel([num2str(max(y_scan),'%.0d'), ' \mum']);
    %     set(ax1,'FontSize',11);
    pbaspect(ax1,[max(x_scan),max(y_scan), 1])
    %     hold off
    [counts,~]=histcounts(scan_image(mask_cell),bin_edges);%,'Normalization','probability');
    [cointsM,~]=histcounts(scan_image(mask_medium),bin_edges);%,'Normalization','probability');
    
    plot(ax3,bin_centers,counts,'-');
    hold(ax3,'on');
    plot(ax3,bin_centers,countsM,'--');
    hold(ax3,'off');
    
end
%% RESIDUAL CODE THAT WAS MOVED TO detect_cells_in_brightfield.m
% %% plot Brightfield
% 
% abs_roi=rect_roi;
% abs_roi(1)=abs_roi(1)+rect_crop(1);
% abs_roi(2)=abs_roi(2)+rect_crop(2);
% 
% rect_roi_full=rect_roi;
% rect_roi_full(1)=rect_roi_full(1)+rect_crop(1);
% rect_roi_full(2)=rect_roi_full(2)+rect_crop(2);
% 
% hf2=figure(2);
% ax2=axes;
% hf2.Position=[550 230 500 400];
% title(['name: ', name, ' brightfield']);
% % h_bf=imagesc(full_brightfield);
% set(ax2,'ydir','normal');
% shading interp;
% colormap(ax2, 'gray');
% axis tight
% axis equal
% set(ax2,'ylim',[rect_crop(1)-100,rect_crop(1)+rect_crop(3)+200]);
% set(ax2,'xlim',[rect_crop(2)-100,rect_crop(2)+rect_crop(4)+200]);
% 
% 
% hr=rectangle('Position',rect_roi_full);
% hr.EdgeColor=1-[0.1 0.1 0.1];
% hr.LineWidth=1;
% set(ax2,'FontSize',11)
% hold off
% 
% 
% 
% 
% %% Select Brightfield shape

% % 
% % %edge detect cell
% % [~,threshold] = edge(full_brightfield,'sobel');
% % fudgeFactor = 0.8;
% % BWs = edge(full_brightfield,'sobel',threshold * fudgeFactor);
% % se90 = strel('line',5,90);
% % se0 = strel('line',5,0);
% % BWsdil = imdilate(BWs,[se90 se0]);
% % BWdfill = imfill(BWsdil,'holes');
% % seD = strel('diamond',5);
% % BWfinal = imerode(BWdfill,seD);
% % BWfinal = imerode(BWfinal,seD);
% % 
% % % display - this must be here for selection in the next code-paragraph
% % h_bf=imagesc(ax2,BWfinal);
% % set(ax2,'ylim',[rect_crop(1)-100,rect_crop(1)+rect_crop(3)+200]);
% % set(ax2,'xlim',[rect_crop(2)-100,rect_crop(2)+rect_crop(4)+200]);
% % colormap(ax2,'gray');
% % 
% % %select by hand only the region of the cell
% % h=imfreehand(ax2);
% % bw_bf=createMask(h,h_bf);
% % bw_bf=logical(bw_bf.*BWfinal);
% % 
% % %display again
% % h_bf=imagesc(ax2,bw_bf);
% % set(ax2,'ylim',[rect_crop(1)-100,rect_crop(1)+rect_crop(3)+200]);
% % set(ax2,'xlim',[rect_crop(2)-100,rect_crop(2)+rect_crop(4)+200]);
% % colormap(ax2,'gray');
% % waitforbuttonpress;
% % 
% % % % simple smoothing - not used 
% % % windowSize = 10;
% % % kernel = ones(windowSize) / windowSize ^ 2;
% % % bw_bf = conv2(single(bw_bf), kernel, 'same');
% % 
% % % %display again
% % % h_bf=imagesc(ax2,bw_bf);
% % % set(ax2,'ylim',[rect_crop(1)-100,rect_crop(1)+rect_crop(3)+200]);
% % % set(ax2,'xlim',[rect_crop(2)-100,rect_crop(2)+rect_crop(4)+200]);
% % % colormap(ax2,'gray');
% % % waitforbuttonpress;
% % 
% % %last line for simple thrsholding - do not use.
% % % bw_bf = bw_bf > 0.5; % Rethreshold
% % delete(h);
% % 
% % 
% % 
% % Bbf=bwboundaries(bw_bf);
% % if length(Bbf)>1
% %     len=zeros(length(Bbf),1);
% %     for k=1:length(Bbf)
% %         len(k)=size(Bbf{k},1);
% %     end
% %     [~,idx_b]=max(len);
% %     Bbf=Bbf(idx_b);
% % end
% % bbf=Bbf{1};
% % 
% % windowWidth = 35;
% % polynomialOrder = 2;
% % smoothX = sgolayfilt(bbf(:,2), polynomialOrder, windowWidth);
% % smoothY = sgolayfilt(bbf(:,1), polynomialOrder, windowWidth);
% % 
% % bbf=[smoothY,smoothX];
% % 
% % imagesc(ax2,full_brightfield);
% % set(ax2,'ylim',[rect_crop(1)-100,rect_crop(1)+rect_crop(3)+200]);
% % set(ax2,'xlim',[rect_crop(2)-100,rect_crop(2)+rect_crop(4)+200]);
% % colormap(ax2,'gray');
% % hold(ax2,'on');
% % plot(bbf(:,2),bbf(:,1),'r-');
% % 
% % hold(ax2,'off');
% 
% hf={hf1,hf2};
% ax={ax1,ax2};
close(hf1);
close(hf3);
end

