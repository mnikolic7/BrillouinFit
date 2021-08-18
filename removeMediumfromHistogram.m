function [yRemoved,average,med_average,bin_centers2] = removeMediumfromHistogram(bin_centers,yhist)
%This function takes in the histogram of the flow/suspended Brillouin shift
%data and takes the first 9 bins as the location where the medium peak is.
%It fits one gaussian to the histogram, and subtracts it from the yhist.
%The output is the remaining histogram of the shifts. it also returns the
%average value of the remaining histogram.
%   author: mnikolic@umd.edu
%%
%make bin_centers a row vector
if ~isrow(bin_centers)
    bin_centers=transpose(bin_centers);
end
%make yhist a column vector
if isrow(yhist)
    yhist=transpose(yhist);
end
%
IDX_M=9;
gaussEqn = 'a*exp(-((x-b)/c)^2)';
startPoints = [0.3 6.06 0.01];
bin_centers_short=bin_centers(1:IDX_M);
yhist_medium=yhist(1:IDX_M);

f = fit(bin_centers_short',yhist_medium ,gaussEqn,'Start', startPoints);
med_average=f.b;
yFit=f(bin_centers);
% yFit=transpose(yFit);
% plot(bin_centers_short,yhist_medium);
% hold on
% plot(f);
% hold off
yRemoved=yhist-yFit;

bin_centers2=bin_centers(IDX_M+1:end);
yRemoved2=yRemoved(IDX_M+1:end);
yRemoved2(yRemoved2<0)=0;
yRemoved2=yRemoved2./sum(yRemoved2); 
% disp 'Returned vector is normalized - in units of probability';
average=bin_centers2*yRemoved2;
% disp(['average shift = ', num2str(average)]);
yRemoved=yRemoved2;
end

