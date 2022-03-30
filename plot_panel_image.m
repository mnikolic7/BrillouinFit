% load('c3.mat');
% load('c3_cells.mat');
% AR=nan(length(cells),1);
% m=nan(length(cells),1);
% circularity=nan(length(cells),1);
BW=~isnan(scan_image2);
BW=bwareafilt(BW,1);
stats=regionprops(BW,'Centroid');
CENTROID = stats.Centroid;
XC=CENTROID(1);
YC=CENTROID(2);
B_large=bwboundaries(BW);
b_large=B_large{1};
%

panel_image=mean(medium)*ones(size(scan_image))*nan;
for k=1:length(cells)
    bw=masks{k};
%     bw=
    panel_image(bw==1)=nanmean(cells{k});
%     stats = regionprops(bw,'Centroid','Area','Perimeter','MajorAxisLength','MinorAxisLength');
    stats=regionprops(bw,'Centroid');
    centroid = stats.Centroid;
    xc(k)=centroid(1);
    yc(k)=centroid(2);
%     circularity(k)=(stats.Area*4*pi)/(stats.Perimeter).^2;
%     AR(k)=stats.MajorAxisLength./stats.MinorAxisLength;
    m(k)=nanmean(cells{k});
    d_center(k)=sqrt((xc(k)-XC)^2+(yc(k)-YC)^2);
    
    db=sqrt((b_large(:,2)-xc(k)).^2+(b_large(:,1)-yc(k)).^2);
    db=min(db);
    d_boundary(k)=db;
end
idx=8;
mean_cells(1:length(m),idx)=m';
D_center(1:length(d_center),idx)=d_center';
D_boundary(1:length(d_boundary),idx)=d_boundary';
clear cells
clear cells_shift
clear masks
clear d_center
clear d_boundary
clear m
%%
%
hf=figure(2);
hf.Position=[50 50 900 550];
ax1=axes();
% ax1=subplot(1,2,1);

% pcolor(x_scan,y_scan,scan_image)
pcolor(ax1,x_scan,y_scan,panel_image);
set(gca,'ydir','reverse');
hold on
% b=B{1};
% plot(x_scan(b(:,2)),y_scan(b(:,1)),'-','color',color_darkred);
for k=1:length(m)
    text(x_scan(round(xc(k))),y_scan(round(yc(k))),num2str(m(k),'%.3f'),'HorizontalAlignment','Center');
end
axis equal
axis tight
shading flat
colormap(ax1, 'jet');
colorbar();
caxis([cmap_min_shift cmap_max_shift]);
axis equal
set(gca,'ytick',[],'xtick',[]);
xlabel([num2str(x_distance,'%.0d'), ' \mum']);
ylabel([num2str(y_distance,'%.0d'), ' \mum']);
title(name);
% title(['Brillouin: z = ', num2str(z_pos-6034.7), '\mum']);
% plot(all_boundary_y2,all_boundary_x2,'--','color',color_lightblue,'linewidth',0.5);
set(gca,'FontSize',11);
%%
saveas(hf,[save_dir,'figures\',name,'_cells.png']);
close(hf)
clear;