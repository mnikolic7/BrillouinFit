function [mask_cell] = BrillouinThreshold(scan_image,mask_cell,mask_medium)
%This function takes in the image (scan_image) of the Brillouin scan, and a
%bw image of the same size that indicates a hand selection of a region away
%from cell (in the medium). It find the mean and std of the shift in the
%medium region and calculates the 3 sigma threshold (99.9978%) for shifts in the
%cell image. Then it returns a mask image that indicates where the cell is.
%It also takes in teh mask_cell image which is a region beyond which
%   author: mnikolic@umd.edu
medium=scan_image(mask_medium);
medium=medium(:);
mean_medium=mean(medium);
std_medium=std(medium);
thresh=mean_medium+3*std_medium;

bw=scan_image>thresh;
mask_cell=bw.*mask_cell;
mask_cell=logical(mask_cell);
end

