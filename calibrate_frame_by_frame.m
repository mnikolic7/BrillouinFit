function [GHz_per_pixel,FSR,get_data_info]=calibrate_frame_by_frame(start_path,varargin)
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
shift_GHz_m=5.57*532/660; %GHz converted to 660 nm from 532 nm laser

%% Open the data files and ask the user to select the peaks. 
select='guess';

% if nargin>=4
%     select=varargin{1}; 
%     y_selected=varargin{2};
%     x_range=varargin{3};
%     rect_bak=varargin{4};
% end
% disp(select)
%% for methanol
[x_m,data2]=get_data_from_file([pathname2, filename2],select);
y_m=mean(data2,2);
%% for water
[x_w,data1,y_selected, x_range, rect_bak]=get_data_from_file([pathname1, filename1],select);
y_w=mean(data1,2);
% %% interpolate
% xf_m=linspace(x_m(1),x_m(end),10*length(x_m));
% xf_w=linspace(x_w(1),x_w(end),10*length(x_w));
% y_m=interp1(x_m,y_m,xf_m,'pchip');
% y_w=interp1(x_w,y_w,xf_w,'pchip');

%%
%save the selected location of the peaks in the water images for optional further processing.
get_data_info=struct('y_selected',y_selected, 'x_range',x_range, 'rect_bak',rect_bak);

%%%%%%%% calculate the GHz_per_pixel %%%%%%%%
hf1=figure(2);

[peak_positions_w,param0]=fitPeaks(x_w,data1(:,1),2, 'guess', [],'displayon');
seps_w(1)=abs(peak_positions_w(2)-peak_positions_w(1));
for k=2:size(data1,2)
    [peak_positions_w,param0]=fitPeaks(x_w,data1(:,k),2, 'provide', param0,'displayoff');
    params_w(:,k)=param0;
%     title('Water fit preview, click/press any key to continue');
    % waitforbuttonpress;
%     plot(x_w,data1(:,k),'.-');
%     pause(0.001);
    seps_w(k)=abs(peak_positions_w(2)-peak_positions_w(1));
end

[peak_positions_m,param0]=fitPeaks(x_m,data2(:,1),2, 'guess', [],'displayoff');
seps_m(1)=abs(peak_positions_m(2)-peak_positions_m(1));
for k=2:size(data2,2)
    [peak_positions_m,param0]=fitPeaks(x_m,data2(:,k),2, 'provide', param0,'displayoff');
    params_m(:,k)=param0;
%     title('methanol fit preview, click/press any key to continue');
    % waitforbuttonpress;
%     plot(x_m,data2(:,k),'.-');
%     pause(0.001);
    seps_m(k)=abs(peak_positions_m(2)-peak_positions_m(1));
end
% watiforbuttonpress;


%%
plot(params_w(1,:)-mean(params_w(1,:)))
hold on
plot(params_w(4,:)-mean(params_w(4,:)))
plot(seps_w-mean(seps_w),'k-');
hold off
waitforbuttonpress;
close(hf1);
%%
sep_w=mean(seps_w);
sep_m=mean(seps_m);
%%
sep_w=abs(peak_positions_w(2)-peak_positions_w(1));
sep_m=abs(peak_positions_m(2)-peak_positions_m(1));

GHz_per_pixel=2*(shift_GHz_w-shift_GHz_m)./(sep_m-sep_w);
FSR=GHz_per_pixel*sep_m+2*shift_GHz_m;
end
