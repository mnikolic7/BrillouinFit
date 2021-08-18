function [h,p] = testNormality(y,varargin)
%This function tests whether the data in vector y is normally distributed
%using kolmogorov smirnov test.
%   author: mnikolic@umd.edu

plotOn=0;
if nargin>1
    plotOn=varargin{1};
end

y=(y-nanmean(y))./nanstd(y);

[h,p]=kstest(y);

if h==0
    disp(['distribution is normal and p value is ',num2str(p),'. P value is the probability that the opposite conclusion (not normal distribution) would be due to chance.']);
else
    disp(['distribution is not normal and there is only ',num2str(p),' probability that this conclusion is due to chance']);
end
if plotOn
    cdfplot(y);
    hold on
    x_values = linspace(min(y),max(y));
    plot(x_values,normcdf(x_values,0,1),'r-')
    legend('Empirical CDF','Standard Normal CDF','Location','best')
    hold off
end
end

