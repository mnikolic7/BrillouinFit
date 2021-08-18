function [radial_profile,mean_radial_profile,r,th] = get_radial_profile(vData)
%This function takes in the data structure outputted by the
%calculate_cell_metrics.m and extracts for each entry (each cell) the
%radial profile of the image. Then it packages them all as a matrix (200
%rows for each radius, by N where n is the number of cells). it also
%returns the mean radial profile (across cells).
%   author: mnikolic@umd.edu

r=vData(1).r_scan;
th=vData(1).th_scan;
N=length(vData);
r=r./max(r); %normalize this to 1. 
radial_profile=nan(200,N);
% figure(1)
% hold on
for k=1:N
%     plot(r,ddd(k).radial_profile,'-');
    radial_profile(:,k)=vData(k).radial_profile;
end
% hold off
%
mean_radial_profile=nanmean(radial_profile,2);

end

