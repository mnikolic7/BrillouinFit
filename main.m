clear;
clc;
%%
warning('off','all');

%%
data_path=('data_path\');

save_dir='save_dir\';
cell_data=struct('name',[],'scan_dir',[],'scan_image',[],...
    'cell',[],'medium',[],'rect_crop',[], 'brightfield', [],...
    'x_scan',[],'y_scan',[],...
    'mask_cell',[],'mask_medium',[],'GHz_per_pixel',[],'FSR',[],...
    'bottom',[],'xy_height',[], 'params',[],'gel',[]);

mkdir(save_dir);
mkdir(save_dir,'figures');
% mkdir(save_dir,'figures_full');
%% calibrate
[GHz_per_pixel, FSR]=calibrate(data_path);

%% get scan file
name='c1';

    


cmap_min_shift=6.0;
cmap_max_shift=6.6;
%
% 
% [filename, pathname] = uigetfile('*.*','Please select the scan data file',...
%     data_path);
% frames=getFrames([pathname,'/',filename]);
acqdatafile=dir([data_path,name,'\*AcqData.tdms']);
pathname=[acqdatafile.folder,'\'];
filename=acqdatafile.name;


[ x, data, y_selected, x_range, rect_bak ]=get_data_from_file([acqdatafile.folder,'\',acqdatafile.name],'manual');
% [ ~, data_refl]=get_data_from_file([acqdatafile.folder,'\',acqdatafile.name],'manual');
%
% R=mean(data_refl,1);

bak=min(data(:));
data=data-bak;
%%%%

% This is the long part: fit the scan
x_fine=linspace(x(1),x(end),200);
axMAX=1.1*max(data(:));
Nframes=size(data,2);
pp1=zeros(Nframes,1);
pp2=zeros(Nframes,1);
params=zeros(6,Nframes);

[peak_positions,param0]=fitPeaks(x,data(:,1),2, 'guess', [],'displayoff');
params(:,1)=param0;
pp1(1)=peak_positions(1);
pp2(1)=peak_positions(2);
%
plotOn=0;
if plotOn
f=figure(2);
% f.Position=[128 512 800 400];
ax=axes;
hold on
plot(ax,x,data(:,1),'bo');
plot(ax,x,curve(param0,x),'r.-');
hold off
end
fprintf(1,['Total frames: ',num2str(Nframes, '%d'), '. Current frame:  000001']);
%
frames_to_fit=2:Nframes;
for n=frames_to_fit
    [peak_positions,param,fval]=fitPeaks(x, data(:,n), 2, 'provide',param0,'displayoff');
%     param0=param;
    pp1(n)=peak_positions(1);
    pp2(n)=peak_positions(2);
    params(:,n)=param;
    if plotOn
    plot(ax,x,data(:,n),'bo');
    hold on
    plot(ax,x,curve(param,x),'r.'); 
    plot(ax,x_fine,curve(param,x_fine),'r-');
    title(ax, [num2str(n),' done out of ', num2str(frames_to_fit(end))]);
    axis([x(1) x(end), 0, axMAX]);
    hold off
    end
    fprintf(1,'\b\b\b\b\b\b%.6d',n); 
    pause(0.0001);
end
fprintf('\n done! \n');

% calculate shifts from fits (just rerun these couple lines with different
% calibrations to compare if needed)

%
seps=abs(pp2-pp1);
seps=seps*GHz_per_pixel;
shifts=sep2shift(seps,FSR); 

widths=GHz_per_pixel*0.5*(params(2,:)+params(5,:));
intensities=0.5*(params(3,:)+params(6,:));
%
%
scan_parameters=readAcqDataInfo([pathname, filename]);
%
z_pos=scan_parameters{'ZPos',1};
scan_rows=scan_parameters{'XCount',1};
scan_columns=scan_parameters{'YCount',1};
scan_columns_z=scan_parameters{'ZCount',1};

x_distance=scan_parameters{'XDis',1};
y_distance=scan_parameters{'YDis',1};
z_distance=scan_parameters{'ZDis',1};


if y_distance<=0
    scan_columns=scan_columns_z;
    y_distance=z_distance;
    scan_type='xz'
elseif x_distance<=0
    scan_rows=scan_columns;
    x_distance=y_distance;
    scan_columns=scan_columns_z;
    y_distance=z_distance;
    scan_type='yz'
else 
    scan_type='xy'
end

bottom=z_pos;

startX=scan_parameters{'Origin Position X',1};
startY=scan_parameters{'Origin Position Y',1};
endX=scan_parameters{'End Position X',1};
endY=scan_parameters{'End Position Y',1};
%


CMOS_umpp=0.1125; %mircons per pixel for 60x objective

%
CMOSroi=[startX, endY, abs(startX-endX), abs(endY-startY)];
%
xmicronpp=x_distance/scan_rows; %x dimension microns per pixel
ymicronpp=y_distance/scan_columns; %y dimension microns per pixel

x_scan=xmicronpp:xmicronpp:scan_rows*xmicronpp;
y_scan=ymicronpp:ymicronpp:scan_columns*ymicronpp;
% %% don't need this section unless we had stopped the scan 
% if length(shifts)<scan_rows*scan_columns
%     shifts(length(shifts)+1:scan_rows*scan_columns)=0;
% end
%
%
files=dir([pathname,'*.tdms']);
for k=1:length(files)
    temp_fname=files(k).name;

    sub_str_idx=strfind(temp_fname,'Green');
    if ~isempty(sub_str_idx)
        bf_fname=temp_fname;
    end
end
%
I_bf=get_image_from_tdms([pathname,bf_fname]);
I_bf=flipud(rot90(I_bf));



xy_height=bottom;
disp(['xy height = ', num2str(xy_height,'%.0d')]);


% I_f=get_image_from_tdms([pathname,f_fname]);
%
margin_size=150;
crop_roi=[CMOSroi(1)-margin_size CMOSroi(2)-margin_size CMOSroi(3)+2*margin_size CMOSroi(4)+2*margin_size];
rect_roi=[margin_size+1 margin_size+1 CMOSroi(3) CMOSroi(4)];
% imagesc(I_bf)
% rect_roi=getrect();

if strcmp(scan_type,'xz')
    rect_roi(4)=1;
elseif strcmp(scan_type,'yz')
    rect_roi(3)=1;
end
%
Ib=imcrop(I_bf,crop_roi);
% Iv=imcrop(I_v,crop_roi);
% If=imcrop(I_f,crop_roi);
% hfbf=figure(6);
% hfbf.Position=[30 30 200 200];
% imshow(Ib,[]);
% title('Please select the cell shape');
% h=imfreehand;
% bw_bf=createMask(h);




%
hf=figure(5);
% get the image of the scan.
%
scan_image=reshape(shifts, scan_rows, []); %scan_image is the result. Just plot this to see the result
scan_image=rot90(scan_image);

scan_image_w=reshape(widths, scan_rows, []); %scan_image is the result. Just plot this to see the result
scan_image_w=rot90(scan_image_w);

scan_image_I=reshape(intensities, scan_rows, []); %scan_image is the result. Just plot this to see the result
scan_image_I=rot90(scan_image_I);
%

% scan_image_r=reshape(refl, scan_rows, []); %scan_image is the result. Just plot this to see the result
% scan_image_r=rot90(scan_image_r);

imagesc(x_scan,y_scan,scan_image);
% pcolor(x_scan,y_scan,scan_image);
%
shading flat
colorbar;
colormap jet
caxis([cmap_min_shift cmap_max_shift]);

title('please select the cell');
h=imfreehand;
bw=createMask(h);
cellData=scan_image(bw);

title('please select the medium');
h=imfreehand;
bw_medium=createMask(h);
medium=scan_image(bw_medium);
% 
scan_image2=scan_image;
scan_image2(~bw)=nan;
thresh=mean(medium)+3.5*std(medium);
scan_image2(scan_image2<thresh)=nan;
%

%
cellData2=cellData;
cellData2(cellData2<thresh)=[];




waitforbuttonpress;
close(hf);



%
% ignore
%
% color_blue=[0, 0.4470, 0.7410];
% color_red=[0.8500, 0.3250, 0.0980];
% color_orange=[0.9290, 0.6940, 0.1250];
% color_purple=[0.4940, 0.1840, 0.5560];  
% color_green=[0.4660, 0.6740, 0.1880];
color_lightblue=[0.3010, 0.7450, 0.9330];
color_darkred=[0.6350, 0.0780, 0.1840];
%
%
bw2=~isnan(scan_image2);
B=bwboundaries(bw2);


if length(B)>1
    len=zeros(length(B),1);
    for k=1:length(B)
        len(k)=size(B{k},1);
    end
    [~,idx_b]=max(len);
    B=B(idx_b);
end


% optional for selecting cells in spheroids
% [LabelMatrix, average_shifts, dist_to_centroid, dist_to_boundary] = select_cells_in_nodules(scan_image,xmicronpp);
        


hf=figure(2);
hf.Position=[50 50 900 550];
ax1=subplot(2,2,1);
%
% hf=figure(3);
% ax1=axes();
% pcolor(x_scan,y_scan,scan_image2)

imagesc(x_scan,y_scan,scan_image);
% set(gca,'ydir','normal');
hold on
b=B{1};
plot(x_scan(b(:,2)),y_scan(b(:,1)),'-','color',color_darkred);
axis equal
axis tight
shading flat
colormap(ax1, 'jet');
colorbar();
caxis([cmap_min_shift cmap_max_shift]); 
% axis equal

set(gca,'ytick',[],'xtick',[]);
xlabel([num2str(x_distance,'%.0d'), ' \mum']);
ylabel([num2str(y_distance,'%.0d'), ' \mum']);
title(['Avg Brillouin shift = ', num2str(nanmean(cellData2),'%.4f')]);
% title(['Brillouin: z = ', num2str(z_pos-6034.7), '\mum']);
% plot(all_boundary_y2,all_boundary_x2,'--','color',color_lightblue,'linewidth',0.5);
% set(gca,'FontSize',11);
pbaspect(ax1,[x_distance,y_distance,1]);
hold off
%
ax2=subplot(2,2,2);
% pcolor(x_scan,y_scan,scan_image)
imagesc(x_scan,y_scan,scan_image_w);
% set(gca,'ydir','normal');
hold on
b=B{1};
% plot(x_scan(b(:,2)),y_scan(b(:,1)),'-','color',color_darkred);

shading flat
colormap(ax2, 'jet');
colorbar();
caxis([0.5 1.5]); 

set(gca,'ytick',[],'xtick',[]);
xlabel([num2str(x_distance,'%.0d'), ' \mum']);
ylabel([num2str(y_distance,'%.0d'), ' \mum']);
title(['Brillouin width']);
% title(['Brillouin: z = ', num2str(z_pos-6034.7), '\mum']);
% plot(all_boundary_y2,all_boundary_x2,'--','color',color_lightblue,'linewidth',0.5);
set(gca,'FontSize',11);
pbaspect(ax2,[x_distance,y_distance,1]);
hold off

ax3=subplot(2,2,3);
imshow(Ib,[]);
% set(gca,'ydir','normal');
shading interp;
hold on
% plot(b_bf(:,2),b_bf(:,1),'-','color',color_lightblue);
colormap(ax3,'gray');
% plot(20+[0 10/CMOSumpp_x]*5 ,[40 40], 'k-','linewidth',6);
axis tight
axis equal
title('Brightfield');
set(gca,'ytick',[],'xtick',[]);

hr=rectangle('Position',rect_roi);
hr.EdgeColor=1-[0.1 0.1 0.1];
hr.LineWidth=1;
set(gca,'FontSize',11)
hold off

ax4=subplot(2,2,4);
% pcolor(x_scan,y_scan,scan_image)
imagesc(x_scan,y_scan,scan_image_I);
% set(gca,'ydir','normal');
hold on
b=B{1};
% plot(x_scan(b(:,2)),y_scan(b(:,1)),'-','color',color_darkred);

shading flat
colormap(ax2, 'jet');
colorbar();
% caxis([200 10e5]); 

set(gca,'ytick',[],'xtick',[]);
xlabel([num2str(x_distance,'%.0d'), ' \mum']);
ylabel([num2str(y_distance,'%.0d'), ' \mum']);

title('Brillouin intensity');
set(gca,'ytick',[],'xtick',[]);
set(gca,'FontSize',11)
pbaspect(ax4,[x_distance,y_distance,1]);
hold off


%
saveas(hf,[save_dir,'figures\',name,'.png']);
%
waitforbuttonpress;
close(hf);
% 
%
cell_data.name=name;
cell_data.scan_dir=pathname;
cell_data.scan_image=scan_image;
cell_data.cell=cellData2;
cell_data.params=params;
cell_data.medium=scan_image(bw_medium);
% cell_data.gel=scan_image(bw_gel);
%
cell_data.brightfield=Ib; 
cell_data.rect_roi=rect_roi;
cell_data.rect_crop=crop_roi;
cell_data.full_brightfield=I_bf;
% cell_data.mask_bf=bw_bf;

cell_data.x_scan=x_scan;
cell_data.y_scan=y_scan;
cell_data.mask_cell=bw2;
cell_data.mask_medium=bw_medium;
cell_data.GHz_per_pixel=GHz_per_pixel;
cell_data.FSR=FSR;
% cell_data.LabelMatrix=LabelMatrix;
% cell_data.cell_shifts=average_shifts;
% cell_data.cell_distToCentroid=dist_to_centroid;
% cell_data.cell_distToBoundary=dist_to_boundary;

cell_data.bottom=bottom;
cell_data.xy_height=xy_height;
%
% str = input('Type the cell filename for saving:   ','s');
save([save_dir,name,'.mat']);
% saveas(hf,[save_dir,'figures\',name,'.png']);