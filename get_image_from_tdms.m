function [data] = get_image_from_tdms( fname )
% This function opens the tdms file given by the directory address fname and returns the 
% single image stored in it as a 2D array. It will also open a tif file and in that case it
% will return all of the frames as a 3D TIFFStack array. 
%
% author: mnikolic@umd.edu
%   dependencies: convertTDMS, TIFFStack

[path,name,ext] = fileparts(fname);
% read tif
if (strcmp(ext, '.tif') || strcmp(ext, '.tiff')) %if the file is a tif file
    frames=TIFFStack(fname); %Use the TIFFStack library (included) to load the file
    [Nrows,Ncolumns,Nframes]=size(frames);
% read tdms
elseif (strcmp(ext, '.tdms')) %if the file is tdms
    %%
    ConvertedData=convertTDMS(0,fname); %load the file using the ConvertTDMS library (included)
    d=ConvertedData.Data.MeasuredData; %extract the data
    %%
    N=length(d); %number of the data columns: vertical (maybe) size of the Brillouin (EMCCD) image plus two empty rows at the beginning (Don't know why).
    Nrows=N-2; %so, this is the number of rows in a single EMCCD image.
    Ncolumns=2560; %Number of horizontal (maybe) size of the Brillouin (EMCCD) image should be the same as the vertical one (number of columns in the actual EMCCD image).
    Nframes=1; %number of frames is the length of each data column, divided by the horizontal size (they are just concatenated together into one data column, again, I don't know why.)
    %%
    frames=zeros(Nrows,Ncolumns,1); %initialize a 3D matrix to save the data in it.
    for n=3:N
        column=d(n).Data; %read the whole data column
        column=reshape(column, [Ncolumns, length(column)./Ncolumns]); %reshape it to make it into a list of rows - one row of horizontal_size for each frame (Nframes total).
        frames(n,:,:)=column; %assign the list of horizontal rows for the vertical position n, to corresponding frames. Think of assigning the horizontal layer of a cube, while the vertical (front) face is the image.
    end
    
else 
    disp 'Currently, this program only supports .tif or .tdms files';
    return
end
data=frames;
end

