function [wc_mean] = sus_whole_cell_mean(list_of_shifts)
%This function takes in the vector of all shifts measured in suspended (or
%flow) experiment, and creates a KDE (uncomment section below for histogram
%option) distribution of those shifts. It removes the medium by fitting the
%first few bins with a Gaussian, and then calculates the average of the
%rest of the distribution.
%Dependencies: create_histogram.m and removeMediumfromHistogram.m. This
%function is basically a wrapper for those two operations.
%   author: mnikolic@umd.edu

%make sure input is a column vector (bootstrap operates on rows)
if isrow(list_of_shifts)
    list_of_shifts=transpose(list_of_shifts);
end
[~,yKDE]=create_histogram(list_of_shifts);
bin_width=0.01; %GHz
bin_centers=6:bin_width:7;
% bin_edges=[bin_centers-bin_width/2 bin_centers(end)+bin_width/2];
[~,wc_mean]=removeMediumfromHistogram(bin_centers,yKDE);

end

