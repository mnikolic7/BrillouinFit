function [GHz_per_pixel,FSR,get_data_info]=calibrate(start_path,varargin)
% This function takes the string directory addres as the input. 
% It prompts the user to select the scan files for methanol and water.
% It calculates the the GHz per pixel value and FSR.
% The third output is the info about the selected lines in the water images.
% so that you can calibrate Brilloin measurements. 
% It averages all the frames in the files that are selected.
% author: mnikolic@umd.edu

%ask the user to select the files.
disp('Please select the water data file');
[filename1, pathname1] = uigetfile('*.*','Please select the water data file',...
    start_path);
disp('Please select the methanol data file now.');
[filename2, pathname2] = uigetfile('*.*','Please select the methanol data file',...
    start_path);

%Literature values for the water and methanol shifts.
shift_GHz_w=7.46*532/660; %GHz converted to 660 nm from 532 nm laser
% shift_GHz_m=5.61*532/660; %GHz converted to 660 nm from 532 nm laser
shift_GHz_m=5.57*532/660; %GHz converted to 660 nm from 532 nm laser

%% Open the data files and ask the user to select the peaks. 
select='manual';

% if nargin>=4
%     select=varargin{1}; 
%     y_selected=varargin{2};
%     x_range=varargin{3};
%     rect_bak=varargin{4};
% end
% disp(select)
%% for methanol
[x_m,data2,y_selected, x_range, rect_bak]=get_data_from_file([pathname2, filename2],select);
y_m=mean(data2,2);
%% for water
[x_w,data1,y_selected, x_range, rect_bak]=get_data_from_file([pathname1, filename1],'provide',y_selected, x_range, rect_bak);
y_w=mean(data1,2);

%%
%save the selected location of the peaks in the water images for optional further processing.
get_data_info=struct('y_selected',y_selected, 'x_range',x_range, 'rect_bak',rect_bak);

%%%%%%%% calculate the GHz_per_pixel %%%%%%%%
hf1=figure(2);
[peak_positions_w,param0]=fitPeaks(x_w,y_w,2, 'guess', [],'displayon');
title('Water fit preview, click/press any key to continue');
waitforbuttonpress;
[peak_positions_m,param0]=fitPeaks(x_m,y_m,2, 'guess', [],'displayon');
title('Methanol fit preview, click/press any key to continue');
waitforbuttonpress;
close(hf1);

%%
sep_w=abs(peak_positions_w(2)-peak_positions_w(1));
sep_m=abs(peak_positions_m(2)-peak_positions_m(1));

GHz_per_pixel=2*(shift_GHz_w-shift_GHz_m)./(sep_m-sep_w);
FSR=GHz_per_pixel*sep_m+2*shift_GHz_m;
% FSR=GHz_per_pixel*sep_w+2*shift_GHz_w;
end
