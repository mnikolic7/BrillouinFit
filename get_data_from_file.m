function [ x, data, y_selected, x_range, rect_bak ] = get_data_from_file( fname, varargin )
%Extracts the brillouin spectrum from tif images and returns an array of
%columns (spectra). One column for each frame.
%
% Two optional arguments: varargin{1} = nothing, manual, or provide
% if manual - will prompt user to find the line where the brillouin
% spectrum lies.
% if provide, must give another argument that is an integer between 1 and
% the vertical number of pixels in the image(s) so that the program will
% pick that line.
%
%dependencies: TIFFStack
% author: mnikolic@umd.edu
%   dependencies: TIFFStack, convertTDMS

% fname
[path,name,ext] = fileparts(fname);
% ext
% read tif
if (strcmp(ext, '.tif') || strcmp(ext, '.tiff')) %if the file is a tif file
    %%
    frames=TIFFStack(fname); %Use the TIFFStack library (included) to load the file
    [Nrows,Ncolumns,Nframes]=size(frames);
    frames=double(frames);
% read tdms
elseif (strcmp(ext, '.tdms')) %if the file is tdms
    %%
    ConvertedData=convertTDMS(0,fname); %load the file using the ConvertTDMS library (included)
    d=ConvertedData.Data.MeasuredData; %extract the data
 
    %%%%%%%%%%%%%%%%%%%%
    %Data is saved like this in the tdms file.

    % d(1).Data=empty;
    % d(2).Data=empty;
    % d(3).Data=[ Nframes*Ncolumns ];
    % d(4).Data=[ Nframes*Ncolumns ];
    % .
    % :
    % d(Nrows).Data=[ Nframes*Ncolumns];
    % 
    % so we need to go for entry of d.Data (each horizontal line in image) and extract the data which is that same row in every frame. 
    % We have to separate the rows from each frame.
    % Save it in the 3-matrix "frames" that is of dimension [Nrows, Ncolumns, Nframes].
    % 
    % Nrows usually = 32
    % Ncolumns usually = 64
    % Nframes usually = varies a lot depending on the scan.

    %%%%%%%%%%%%%%%%%%%%%%
    %%
%     fname=[pathname, '/', filename];
    N=length(d); %number of the data entries: vertical size of the Brillouin (EMCCD) image plus two empty rows at the beginning (Don't know why).
    Nrows=N-2; %so, this is the number of rows in a single EMCCD image.
    Ncolumns=64; %Number of horizontal size of the Brillouin (EMCCD) image (number of columns in the actual EMCCD image).
    Nframes=length(d(3).Data)./Ncolumns; %number of frames is the length of each data column, divided by the horizontal size (they are just concatenated together into one data column, again, I don't know why.)
    %%
    try
    %initialize a 3D matrix to save the data in it.
    frames=zeros(Nrows,Ncolumns,Nframes); %if there is an error here check line 29 and 30 (two lines above):  the width of the EMCCD image in the tdms file is sort of hardwired to be the same as the height, so if it is different (by mistake) you should change it here to whatever number it is in your .tdms files. Don't forget to return Ncolumns to the same value as Nrows because that's the most common case. You don't want to have to worry about this all the time.
    for n=3:N %for each entry of d(n).Data
        row=d(n).Data; %read the whole data row (frames concatenated)
        row=reshape(row, [1, Ncolumns, Nframes]); %reshape it to make it into a list of rows - one row of horizontal_size for each frame (Nframes total).
        frames(n-2,:,:)=row; %assign the list of horizontal rows for the vertical position n, to corresponding frames. Think of assigning the horizontal layer of a cube, while the vertical (front) face is the image.
    end
%     imagesc(frames(:,:,1));
    catch ME
       msg ='Error in reshaping the EMCCD image. Check get_data_from_file.m : lines 53-55. You might need to adapt the numbers to image size.'; 
       causeException = MException('MATLAB:BrillouinFit_v12:EMCCDImageRead',msg);
       ME = addCause(ME,causeException);
       disp(Nrows);
       disp(Ncolumns);
       disp(Nframes);
       rethrow(ME)
    end
elseif (strcmp(ext,'.mat')) %if the file is .mat
    var_name=who('-file',fname);
    load(fname);
    frames=eval(var_name{1});
%     frames=frames_deconv;    
%     [Nrows,Ncolumns,Nframes]=size(frames);
%     vars=who;
%     length(vars)
%     frames=eval(vars{1});
    [Nrows,Ncolumns,Nframes]=size(frames);
else 
    disp 'Currently, this program only supports .tif, .tdms, or .mat files.';
    disp '.mat files must have only one variable that contains the frames of the scan';
    return
end
% frames=frames(1:50,25:75,:);
%%
final_x1=1;
final_x2=Ncolumns;
y_selected=-1;
x_range=-1;

% size(frames)
%% crop the data.
if (nargin>=2 && strcmpi(varargin{1},'manual'))
    %% Continue with the rest. 
    I=frames(:,:,1);

    Good='No';
    while (strcmpi(Good,'No'))
        f = figure(4);
        screensize=get(groot,'Screensize');
        set(f, 'Position', [80 80 screensize(3)/2 screensize(4)-160]); %set the figure position and size
        ax=subplot(2,1,1);
        set(ax,'outerposition',[0 0.5 0.99 0.49]);
        imagesc(ax, I);
        axis equal
        title(ax,'Please select a rectangle that contains the Brillouin signal');
        hold on
        rect=getrect(ax); %G37 R3k7, yo!
        x1=round(rect(1));
        x2=round(rect(3)+rect(1));
        y1=round(rect(2));
        y2=round(rect(2)+rect(4));
        axis([x1-5 x2+5 y1-5 y2+5]);
        rectangle('position',rect);
        title(ax, 'Please click the middle of the Brillouin peaks');
        [~,y_line]=ginput(1);
        y_line=round(y_line);

        title(ax, 'Please select a region with background');
        rect_bak=getrect(ax);
        bak=mean2(I(rect_bak(2):rect_bak(2)+rect_bak(4),rect_bak(1):rect_bak(1)+rect_bak(3)));
        %
        ax2=subplot(2,1,2);
        set(ax2,'outerposition', [0 0 0.99 0.49]);
        profile=mean(I(y_line-2:y_line+2,x1:x2),1);
        plot(ax2, x1:x2, profile)
        title(ax2, 'Horizontal profile selected in the square.');
        
        plot(ax, [x1, x2], [y_line, y_line], 'r-','linewidth', 1);
        plot(ax, [x1, x2], [y_line-2, y_line-2], 'r:','linewidth', 1);
        plot(ax, [x1, x2], [y_line+2, y_line+2], 'r:','linewidth', 1);
        Good = questdlg('Is the selection good?', ...
            'Initial selection good?','Yes','No','No');
        final_x1=x1;
        final_x2=x2;
        y_selected=y_line;
        x_range=[x1, x2];
    end
%     title(ax,'Click to close and proceed');
%     waitforbuttonpress;
    close(f);
    %%
    frames=frames(y_line-2:y_line+2,final_x1:final_x2,:);    
elseif (nargin>=2 && strcmpi(varargin{1},'provide'))
    I=frames(:,:,1);
    y_line=varargin{2};
    x_range=varargin{3};
    
    final_x1=x_range(1);
    final_x2=x_range(2);
    frames=frames(y_line-2:y_line+2,final_x1:final_x2,:);

    rect_bak=varargin{4};
    bak=mean2(I(rect_bak(2):rect_bak(2)+rect_bak(4),rect_bak(1):rect_bak(1)+rect_bak(3)));
else
    disp('Choosing the horizontal line where the Brillouin light is automatically. No Deconvolution'); 
    [~,idx]=max(mean(frames(:,:,1),2)); %crop around the maximum valu
    if Nrows>5
        frames=frames(idx-2:idx+2,:,:);
    end
    bak=200;
    x_range=[1,size(frames,2)];
    y_selected=idx;
    rect_bak=[1 1 1 1];
end

%%
[~,Ncolumns,Nslices]=size(frames);
data=sum(frames,1); %sum along the rows (vertically);
%data=double(data);
%reshape data: each column is the spectrum from one slice, there are
%   Nslices columns
data=reshape(data, Ncolumns, Nslices);
x=1:size(data,1);
x=x';
data=data-bak;

bak=min(data(:));
data=data-bak;
end

