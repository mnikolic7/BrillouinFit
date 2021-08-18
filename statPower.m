function [pwrcurr] = statPower(Data)
%this function takes in a matrix with two columns and uses sampsizepwr
%function to estimate the statistical power of a two sample t test to
%distinguish these two measurements given the number of samples. It does
%not handle nans well yet and if the two conditions have different number of
%samples. It only looks at the first column for N.
%   author: mnikolic@umd.edu
%%
stds=nanstd(Data);
avgs=nanmean(Data);
sigma=mean(stds);
delta=avgs(1)-avgs(2);
N=length(Data(:,1));
%%
% nout = sampsizepwr('t',[avgs(1) stds(1)],avgs(2),0.95)
nn = 1:100;
% pwrout = sampsizepwr('t',[avgs(1) stds(1)],avgs(2),[],nn);
pwr_in=0.95;
nout = sampsizepwr('t2',[avgs(1) sigma],avgs(2),pwr_in);
pwrout = sampsizepwr('t2',[avgs(1) sigma],avgs(2),[],nn);

pwrcurr=sampsizepwr('t2',[avgs(1) sigma], avgs(2),[],N);

figure(1);
plot(nn,pwrout,'b-',nout,pwr_in,'ro')
hold on
plot(N,pwrcurr,'rs');
title({'Power versus Sample Size','given this small sample mean and std',...
    [num2str(nout), ' samples needed for ', num2str(pwr_in) ,' power'],...
    ['Current stat. power = ', num2str(pwrcurr,'%.2f'), ' with N= ', num2str(N)]})
xlabel({'Sample Size','# cells needed'})
ylabel('Statistical Power')
grid on
hold off
end

