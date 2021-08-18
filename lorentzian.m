function [y]=lorentzian(x,x0,wL)
%This function returns the lorentzian profile of width wL, centered at x0,
%evaluated at x
%   wL is equal to gamma of the Lorentzian: twice the FWHM
%
%   Dependencies: none.
%   autor: mnikolic@umd.edu

x=x-x0;
y=1./pi*((wL/2)./((wL/2).^2+x.^2));

return;
end
