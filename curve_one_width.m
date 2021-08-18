function y = curve_one_width(p,x)
% function that is used by fitPeaks. 
% curve that is going to be fiteed to the data.
% author: mnikolic@umd.edu
global Np;
N=Np;
y=zeros(size(x));
wL=p(2);
for k=1:3:N*3
    x0=p(k);
%     wL=p(k+1);
    A=p(k+2);
    y=y+A*lorentzian(x,x0,wL);
end
end