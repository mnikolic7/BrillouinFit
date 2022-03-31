clear;
clc;
%%
load('Oct23_m2_d5_tracked')
load('Oct23_m1_d5.mat');
load('Oct23_m3_d5.mat');
load('Oct23_m4_d5.mat');

%%
vOct23_m2_d5=calculate_cell_metrics(Oct23_m2_d5_tracked);


vOct23_m1_d5=calculate_cell_metrics(Oct23_m1_d5);
vOct23_m3_d5=calculate_cell_metrics(Oct23_m3_d5);
vOct23_m4_d5=calculate_cell_metrics(Oct23_m4_d5);
%%
statPower(D)

%%
shifts=getDataMatrix('mean',vOct23_m1_d5,vOct23_m2_d5,vOct23_m3_d5,vOct23_m4_d5);
linewidths=getDataMatrix('linewidth',vOct23_m1_d5,vOct23_m2_d5,vOct23_m3_d5,vOct23_m4_d5);

%%
groups={'M1 d0','M2 d0','M1 d2','M2 d2','M1 d4','M2 d4','M1 d5','M2 d5'};
cols=[1 6 1 6 1 6 1 6];
plot_bar_scatter(shifts,groups,cols,0);
hold on
plot_bar_scatter_KDE(shifts,groups,cols,0);
ylim([6.2 6.6]);
ylabel('Shift (GHz)');
hold off
%%
%%
groups={'M1 d0','M2 d0','M1 d2','M2 d2','M1 d4','M2 d4','M1 d5','M2 d5'};
cols=[1 6 1 6 1 6 1 6];
plot_bar_scatter(linewidths,groups,cols,1);
hold on
plot_bar_scatter_KDE(linewidths,groups,cols,0);
ylim([0.7 1.5]);
ylabel('Linewidth (GHz)');
hold off
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
stds=nanstd(shifts);
avgs=nanmean(shifts);
sigma=mean(stds);
delta=avgs(1)-avgs(2);
%%
% nout = sampsizepwr('t',[avgs(1) stds(1)],avgs(2),0.95)
nn = 1:100;
% pwrout = sampsizepwr('t',[avgs(1) stds(1)],avgs(2),[],nn);
pwr_in=0.95;
nout = sampsizepwr('t2',[avgs(1) sigma],avgs(2),pwr_in)
pwrout = sampsizepwr('t2',[avgs(1) sigma],avgs(2),[],nn);

figure(1);
plot(nn,pwrout,'b-',nout,pwr_in,'ro')
title({'Power versus Sample Size','given this small sample mean and std',[num2str(nout), ' cells needed for ', num2str(pwr_in) ,' power']})
xlabel({'Sample Size','# cells needed'})
ylabel('Statistical Power')
grid on