function [BW] = detect_cells_in_brightfield(full_brightfield,axlim,rect_roi_full,figure_title)
%DETECT_CELLS_IN_BRIGTFIELD takes in the grayscale image file.
%It detects the edges in the image using the edge detection protocol
%described in the Matlab tutorial https://www.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html
%The function is interactive and allows the user
%to manually edit the selection. 
%The function also takes in figure_title which is just the title of the
%figure in which you will be making the selection. It is useful to display
%the filename there, so that you know which image you are working on.
%
%The function returns: the binarized image containing the selected cells.
%ran on MATLAB R2018b with image analyisis toolbox;
%   author: Milos Nikolic mnikolic@umd.edu mnikolic0706@gmail.com

%% this section is adjusted and taken from: https://www.mathworks.com/help/images/detecting-a-cell-using-image-segmentation.html


BWfinal=edgeDetectCell(normalizeImage( full_brightfield ));

%%

BW=BWfinal;



hf=figure(2);
hf.Position=[400 45 900 1000]; %change this to change position of figure on screen
ax = axes(hf);
ax.Units = 'pixels';
ax.Position = [75 75 750 500]; %change this to change position of image within figure
BFrgb=zeros(size(BWfinal));
BFrgb(:,:,1)=0.75*full_brightfield;
BFrgb(:,:,2)=0.75*full_brightfield;
BFrgb(:,:,3)=0.75*full_brightfield;

MaxI=max(full_brightfield(:));
BFrgb=BFrgb./MaxI;

BWrgb=BFrgb;
BWrgb(:,:,1)=BWrgb(:,:,1)+BWfinal*0.25;

h_bf=imagesc(ax,BWrgb);
% bbf=get_boundary(BW);
hold(ax,'on');
% plot(ax, bbf(:,2),bbf(:,1),'g-');
% hold(ax,'off');
% set(ax,'ydir','normal');
shading interp;
axis tight
axis equal
axis(ax,axlim);

hr=rectangle(ax,'Position',rect_roi_full);
hr.EdgeColor=1-[0.1 0.1 0.1];
hr.LineWidth=1;
% set(ax,'FontSize',11)
hold(ax,'off');

%% Edit selection manually 
%define buttons
c1 = uicontrol(hf);
c1.Position=[110,20,120,30]; %button location
c1.String = 'Select cell region';
c1.Callback = @selectCellRegion;

c2 = uicontrol(hf);
c2.Position=[250,20,120,30]; %button location
c2.String = 'Reset';
c2.Callback = @reset;

c4 = uicontrol(hf);
c4.Position=[390,20,120,30]; %button location
c4.String = 'Add region';
c4.Callback = @addRegionButtonPushed;

c5 = uicontrol(hf);
c5.Position=[530,20,120,30]; %button location
c5.String = 'Close holes+shrink';
c5.Callback = @closeHoles;

c3 = uicontrol(hf);
c3.Position=[670,20,120,30]; %button location
c3.String = 'Finish';
c3.Callback = @endButtonPushed;

%display the second small helper figure
hf2=figure(4);
hf2.Position=[45 45 300 300]; %position of the small figure
ax2=axes(hf2);
imagesc(ax2,BWfinal);
title(ax2,{'Initial guess from edge detection', 'Close only after done with manual selection'});
axis(ax2,axlim);

uiwait(hf2); %wait until the current figure is deleted to update final result 
%(wait until the manual slection is finished).
%return
BW=BWfinal;
return;

%% click and select any regions that contain cells
    function selectCellRegion(src,event)
        title({figure_title,'Draw a region that you want included only',''});
        h=imfreehand;
        BWsel=createMask(h);
        BWfinal(~BWsel)=0;
        BWfinal=smoothMask(BWfinal);
        BW=BWfinal;
        delete(h);
        
        imagesc(ax2,BWfinal);
        title(ax2,{'Selected regions', 'Close only after done with manual selection'});
        axis(ax2,axlim);


        BWrgb=BFrgb;
        BWrgb(:,:,1)=BWrgb(:,:,1)+BWfinal*0.25;
        delete(h_bf)
        
        yl=ylim(ax);
        xl=xlim(ax);
        axlim=[xl, yl];
        %plot from here
        h_bf=imagesc(ax,BWrgb);
        hold(ax,'on');
        hr=rectangle(ax,'Position',rect_roi_full);
        hr.EdgeColor=1-[0.1 0.1 0.1];
        hr.LineWidth=1;

        bbf=get_boundary(BWfinal);
        
        plot(ax, bbf(:,2),bbf(:,1),'g-');
        shading interp;
        axis tight
        axis equal
        hold(ax,'off');

        axis(ax,axlim);

    end

    function reset(src,event)
        BWfinal=edgeDetectCell(normalizeImage(full_brightfield));
        BW=BWfinal;
        
        imagesc(ax2,BWfinal);
        title(ax2,{'Selected regions', 'Close only after done with manual selection'});
        axis(ax2,axlim);
        
        BWrgb=BFrgb;
        BWrgb(:,:,1)=BFrgb(:,:,1)+BWfinal*0.25;
        
        delete(h_bf);
        
        yl=ylim(ax);
        xl=xlim(ax);
        axlim=[xl, yl];
        %plot from here
        h_bf=imagesc(ax,BWrgb);
        hold(ax,'on');
        hr=rectangle(ax,'Position',rect_roi_full);
        hr.EdgeColor=1-[0.1 0.1 0.1];
        hr.LineWidth=1;
        shading interp;
        axis tight
        axis equal
        hold(ax,'off');
%         set(ax,'xtick',[],'ytick',[]); 
        axis(ax,axlim);
    end

    function addRegionButtonPushed(src,event)
        h=imfreehand;
        BWsel=createMask(h);
        BWfinal(BWsel)=1;
        delete(h);
        BW=BWfinal;
        

        imagesc(ax2,BWfinal);
        title(ax2,{'Selected regions', 'Close only after done with manual selection'});
        axis(ax2,axlim);
                
        BWrgb=BFrgb;
        BWrgb(:,:,1)=BFrgb(:,:,1)+BWfinal*0.25;
        delete(h_bf);
        
        yl=ylim(ax);
        xl=xlim(ax);
        axlim=[xl, yl];
        %plot from here
        h_bf=imagesc(ax,BWrgb);
        hold(ax,'on');
        hr=rectangle(ax,'Position',rect_roi_full);
        hr.EdgeColor=1-[0.1 0.1 0.1];
        hr.LineWidth=1;

        bbf=get_boundary(BWfinal);
        shading interp;
        axis tight
        axis equal
        plot(ax, bbf(:,2),bbf(:,1),'g-');
        hold(ax,'off');
%         set(ax,'xtick',[],'ytick',[]);
        axis(ax,axlim);
    end

    function closeHoles(src,event)

        BWfinal=imfill(BWfinal,'holes');
        seD = strel('diamond',5);
        BWfinal = imerode(BWfinal,seD);
        BWfinal = imerode(BWfinal,seD);
        BWfinal=smoothMask(BWfinal);
        BW=BWfinal;
        

        imagesc(ax2,BWfinal);
        title(ax2,{'Selected regions', 'Close only after done with manual selection'});
        axis(ax2,axlim);
                
        BWrgb=BFrgb;
        BWrgb(:,:,1)=BFrgb(:,:,1)+BWfinal*0.25;
        delete(h_bf);
        
        yl=ylim(ax);
        xl=xlim(ax);
        axlim=[xl, yl];
        %plot from here
        h_bf=imagesc(ax,BWrgb);
        hold(ax,'on');
        hr=rectangle(ax,'Position',rect_roi_full);
        hr.EdgeColor=1-[0.1 0.1 0.1];
        hr.LineWidth=1;

        bbf=get_boundary(BWfinal);
        shading interp;
        axis tight
        axis equal
        plot(ax, bbf(:,2),bbf(:,1),'g-');
        hold(ax,'off');
%         set(ax,'xtick',[],'ytick',[]);
        axis(ax,axlim);
    end

    function endButtonPushed(src,event)
        BW=BWfinal;

        
        imagesc(ax2,BW);
        title(ax2,{'Final bw mask image', 'You may close this figure now'});
        axis(ax2,axlim);
        close(hf);
    end

    function bbf=get_boundary(mask)
        Bbf=bwboundaries(mask);
        if length(Bbf)>1
            len=zeros(length(Bbf),1);
            for k=1:length(Bbf)
                len(k)=size(Bbf{k},1);
            end
            [~,idx_b]=max(len);
            Bbf=Bbf(idx_b);
        end
        
        bbf=Bbf{1};
        %smooth boundaries
        windowWidth = 11;
        polynomialOrder = 2;
        smoothX = sgolayfilt(bbf(:,2), polynomialOrder, windowWidth);
        smoothY = sgolayfilt(bbf(:,1), polynomialOrder, windowWidth);      
    end

    function mask=boundary2mask(b,I)
        mask=zeros(size(I));
        mask(round(b(:,2)),round(b(:,1)))=1;
        mask=imfill(mask,'holes');
    end

    function mask=smoothMask(mask)
        windowSize = 15;
        kernel = ones(windowSize) / windowSize ^ 2;
        mask = conv2(single(mask), kernel, 'same');
        mask = mask > 0.5; % Rethreshold
    end

    function mask=shrinkMask(mask,Npx)
        for k=1:Npx
            Bbf=bwboundaries(mask);
            bbf=Bbf{1};
            mask(round(bbf(:,2)),round(bbf(:,1)))=0;
        end
    end

    function BWout=edgeDetectCell(brightfield)
        %edge detect cell
        [~,threshold] = edge(brightfield,'sobel');
        fudgeFactor = 0.85;
        BWs = edge(brightfield,'sobel',threshold * fudgeFactor);
        se90 = strel('line',5,90);
        se0 = strel('line',5,0);
        BWsdil = imdilate(BWs,[se90 se0]);
%         BWsdil = imfill(BWsdil,'holes');
        BWout=BWsdil;
%         seD = strel('diamond',5);
%         BWout = imerode(BWsdil,seD);
%         BWout = imerode(BWout,seD);
%         BWout=smoothMask(BWout);
    end
end

