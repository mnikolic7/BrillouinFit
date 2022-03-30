function [LabelMatrix, average_shifts, dist_to_centroid, dist_to_boundary] = select_cells_in_nodules(scan_image,pixel_size)
% This is a helper function that takes in the brillouin image of a
% spheroids and asks you repeatedly to select all the regions of the cells.
% It plots the boundary and the centroid of the whole spheroid, and then allows
% you to manually select the region of the cells until you are satisfied and  
% confirm that you are finished. It returns the label matrix corresponding
% to the each cell selected, positions of the centroid and the boundary,
% and distance from centroid and boundary of each cell.  
% also it is optional to repeat an image, or repeat the image but only
% manually select centroid.
%
% Additionally the function takes in the pixel size so that it can return
% the distances to centroid and boundary in microns instead of pixels.
%
%author: mnikolic@umd.edu

%%
LabelMatrix=zeros(size(scan_image));
Good='Yes';
n=1;
t=6.3;
bw=scan_image>t;
bw=bwareafilt(bw,1);
CC=regionprops(bw,'Centroid');
B=bwboundaries(bw);
centroid=CC.Centroid;


hf=figure(2);
hf.Position=[50 50 900 574];
% hf.Position=[-1749,-145,992,776];
% disp('The figure is on the second monitor, if you are not connected, you won''t see the figure');
ax1=axes();
imagesc(ax1,scan_image);
hold on
b=B{1};
plot(b(:,2),b(:,1),'-','color','k','linewidth',1);
plot(centroid(2),centroid(1),'go','markerfacecolor','g');
axis equal
axis tight
shading(ax1,'flat')
colormap(ax1, 'jet');
colorbar();
caxis([6.0 6.6]);
set(gca,'ytick',[],'xtick',[]);

manual_centroid = questdlg('Select centroid manually?', ...
    'Select centroid manually?','Yes','No','No');
if strcmp(manual_centroid,'Yes')
    [xc,yc]=ginput(1);
    centroid=[yc,xc];
    plot(centroid(2),centroid(1),'ms','markerfacecolor','m');
end
%

while (strcmpi(Good,'Yes'))
    
    title('please select a cell');
    h=imfreehand;
    bw_single=createMask(h);
    bw_single=bw_single & bw;
    LabelMatrix(bw_single)=n;
    
    si=scan_image;
    si(~bw_single)=nan;
    
    B_sc=bwboundaries(bw_single);
    b_sc=B_sc{1};
    plot(ax1,b_sc(:,2),b_sc(:,1),'-','color',[0.7 0.1 0.1],'linewidth',1);
    
    retry = questdlg('Retry?', ...
        'Retry previous cell?','Yes','No','No');
    if strcmp(retry,'Yes')
        n=n-1;
    end
    Good = questdlg('Another?', ...
        'Select another cell?','Yes','No','Yes');
    clear h
    n=n+1;
end
hold off

close(hf);

% calculate some properties to return, but this can also be done later if you also save the label matrix. 
Ncells=n-1;
stats=regionprops(LabelMatrix,scan_image,'Centroid','MeanIntensity');
average_shifts=[stats.MeanIntensity];
dist_to_centroid=nan(Ncells,1);
dist_to_boundary=nan(Ncells,1);

for k=1:Ncells
cell_c=stats(k).Centroid;
    cell_x=cell_c(1);
    cell_y=cell_c(2);
    dist_to_centroid(k)=sqrt((cell_x-centroid(2)).^2+(cell_y - centroid(1)).^2);
    
    
    min_dist=min(sqrt((b(:,2)-cell_x).^2+(b(:,1)-cell_y).^2));
    dist_to_boundary(k)=min_dist;
end

dist_to_boundary=dist_to_boundary*pixel_size;
dist_to_centroid=dist_to_centroid*pixel_size;

% %% check code for the calculation of distance to boundary and distance to
% % centroid.
% 
% hf=figure(1);
% hf.Position=[19.6666666666667,61,1252,571.333333333333];
% ax2=subplot(1,2,2);
% ax1=subplot(1,2,1);
% imagesc(ax1,scan_image)
% colormap jet
% axis equal
% hold on
% plot(centroid(2),centroid(1),'gx');
% plot(ax1,b(:,2),b(:,1),'-','color','k','linewidth',1);
% plot(ax1,b(1,2),b(1,1),'rs','markerfacecolor','r','linewidth',1);
% for k=1:Ncells
%     cell_c=stats(k).Centroid;
%     cell_x=cell_c(1);
%     cell_y=cell_c(2);
%     plot(cell_x,cell_y,'b*');
%     plot([cell_x,centroid(2)],[cell_y centroid(1)],'b-')
%     [~,idx]=min(sqrt((b(:,2)-cell_x).^2+(b(:,1)-cell_y).^2));
%     
%     plot([cell_x,b(idx,2)],[cell_y b(idx,1)],'k-')
%     plot(b(idx,2),b(idx,1),'ms','markerfacecolor','m','markersize',12);
%     plot(ax2,sqrt((b(:,2)-cell_x).^2+(b(:,1)-cell_y).^2))
%     hold(ax2,'on');
%     plot(ax2,idx,sqrt((b(idx,2)-cell_x).^2+(b(idx,1)-cell_y).^2),'ro');
%     hold(ax2,'off');
%     waitforbuttonpress;
% end
% 
% hold off

end
