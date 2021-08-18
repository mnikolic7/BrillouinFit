function [ frames ] = getFrames( fname )
%This function is a part of the get_data_from_file - but it doesn't do
%anything fancy, except that it just gets the RAW EMCCD image and returns
%it as an array of 16bit matrices - a 3D matrix
%   Dependencies: TIFFStack and ConvertTDMS

%% Part 1 - copied from get data from file so that we can obtain all the EMCCD frames 
[path,name,ext] = fileparts(fname);
% read tif
if (strcmp(ext, '.tif') || strcmp(ext, '.tiff')) %if the file is a tif file
    %%
    frames=TIFFStack(fname); %Use the TIFFStack library (included) to load the file
    [Nrows,Ncolumns,Nframes]=size(frames);
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
    % Nrows usually = 10
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
       msg ='Error in reshaping the EMCCD image. Check get_data_from_file.m : lines 43-45. You might need to adapt the numbers to image size.'; 
       causeException = MException('MATLAB:BrillouinFit_v12:EMCCDImageRead',msg);
       ME = addCause(ME,causeException);
       disp(Nrows);
       disp(Ncolumns);
       disp(Nframes);
       rethrow(ME)
    end
elseif (strcmp(ext,'.mat')) %if the file is .mat
    load(fname);
    var=who;
    frames=eval(var{1});
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
end

