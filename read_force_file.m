function [M_ex,M_re] = read_force_file(file)
%This function reads the processed curve file from the JPK AFM analysis
%program that saves the data in a text file. Input must be full path to the
% .txt file containing the measured force displacement and other data from
% the JPK AFM system. It reads the matrix, returns two N x 8 matrices.
% First is the extension part and the second is the retraction part. 8
% columns are: "Vertical Tip Position" "Vertical Deflection" "Height" 
%              "Fast Scanner Height" "Height (measured & smoothed)" 
%              "Height (measured)" "Series Time" "Segment Time"
%   author: mnikolic@umd.edu
M=readmatrix(file);

M=M(:,1:8);
idx=isnan(M(:,1));
IDX_START=find(diff(idx)==-1);
IDX_END=[find(diff(idx)==1),length(idx)];

M_ex=M(IDX_START(1):IDX_END(1),:);
M_re=M(IDX_START(2):IDX_END(2),:);
end

