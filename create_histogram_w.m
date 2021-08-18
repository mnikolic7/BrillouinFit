function [yCounts,yKDE] = create_histogram_w(shift_data)
%CREATE_HISTOGRAM This function takes in a vector of data of Brillouin
%shifts and creates a histogram and returns the values binned in default
%bins, so that you can plot easily after. 
%   author: mnikolic@umd.edu
bin_width=0.01; %GHz
bin_centers=0.5:bin_width:2;
bin_edges=[bin_centers-bin_width/2 bin_centers(end)+bin_width/2];


[yCounts,~]=histcounts(shift_data,bin_edges,'Normalization','probability');
pdKDE = fitdist(shift_data,'kernel','BandWidth',bin_width);
yKDE=pdf(pdKDE, bin_centers)*bin_width;

end

