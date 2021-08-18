function y = curve_providePDMS(p,x)
% function that is used by fitPeaks. 
% curve that is going to be fiteed to the data.
% author: mnikolic@umd.edu

% global Np;
global existing_spectrum_param;

x0_ex1=existing_spectrum_param(1);
x0_ex2=existing_spectrum_param(4);
wL_ex1=existing_spectrum_param(2);
wL_ex2=existing_spectrum_param(5);

y=zeros(size(x)); 
Aex_1=p(7);
Aex_2=p(8);

y=y+Aex_1*lorentzian(x,x0_ex1,wL_ex1);
y=y+Aex_2*lorentzian(x,x0_ex2,wL_ex2);

% fit parameters for cell/medium/water are 
x0_1=p(1);
x0_2=p(4);

wL_1=p(2);
wL_2=p(5);

A_1=p(3);
A_2=p(6);

y=y+A_1*lorentzian(x,x0_1,wL_1);
y=y+A_2*lorentzian(x,x0_2,wL_2);

end