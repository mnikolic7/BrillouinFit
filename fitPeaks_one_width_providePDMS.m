function [ peak_positions,param,fval ] = fitPeaks_one_width_providePDMS( x, y, Npeaks, method, peaks0, existing_spectrum, varargin )
%This function fits Npeaks lorentzian peaks to the data provided in the
%vector y. Make sure x and y have the same dimensions.
% it takes in the x and y data (what you would use to plot the Brillouin spectrum),
% Npeaks is the number of Brillouin peaks (should almost always be 2)
% method is the method in which the initial guess will be selected. 
%        default is 'guess' (could be anything it just has to be a string), for which
%               the program will guess automatically where the peaks are.
%        you can also use 'provide' which will use the parameters for the fit provided 
%               in the next input argument, or
%        'ask' which will ask the user for graphical input for the peaks. 
% peaks0 is the initial guess for the fit parameters, provided that the method is 'provide'
% if you supply any additional input arguments and the last one of them is the string
%   'displayon' in that case this function will also plot the result in the figure 2.
% the output arguments are peak_positions (a Npeaks  by 1 vector of peak positions),
% param (the best fit parameters of the curve, and fval which is the squared error of the fit)
% param is an 3*Npeaks by 1 vector, that contains the center, width, and intensity for each peak 
% in that order: [center1; width1; intenstiy1; center2, width2; intensity2; etc.] 
%   author: mnikolic@umd.edu
%%
Npeaks=2;
method='ask';
peaks0=[];
existing_spectrum=param_PDMS_averaged;
%%
fun=@curve_providePDMS;
intensity_constant=pi;
global Np;
Np=Npeaks;
global existing_spectrum_param;
existing_spectrum_param=existing_spectrum;
%%
plotON=0;
if (nargin>5) && strcmpi(varargin{1},'displayon')
    plotON=1;
end
    
%%
%find peaks for the initial guess of where the centers of the peaks are
[ypeaks, yp_indx]=findpeaks(smooth(y),'NPeaks',Npeaks,'SortStr','descend');
% yp_indx
% figure(1)
% plot(x,y,'.-');
% hold on
% plot(x(yp_indx),ypeaks,'o');
% hold off

if (length(ypeaks)<Npeaks)
%     disp 'error! is here. Found only one peak in fitPeaks.m'
%     method='ask';
    ypeaks=repmat(max(y),[1,Npeaks]);
    yp_indx=round(linspace(1,length(x),Npeaks));
end

param0=zeros(1,3*Npeaks);
for n=1:Npeaks
    idx_x0=(n-1)*3+1;
    idx_wL=(n-1)*3+2;
    idx_A=(n-1)*3+3;
    
    param0(idx_x0)=x(yp_indx(n));
    param0(idx_wL)=3;
    param0(idx_A)=ypeaks(n)*intensity_constant;
end
param0_extended=param0;
param0_extended(length(param0)+1)=existing_spectrum(3);
param0_extended(length(param0)+2)=existing_spectrum(6);

%%%%%%%%%%%%%%%%
if strcmpi(method, 'provide')
    param0=peaks0;
    param0_extended=param0;
    param0_extended(length(param0)+1)=existing_spectrum(3);
    param0_extended(length(param0)+2)=existing_spectrum(6);
%%%%%%%%%%%%%%%%
elseif strcmpi(method,'ask')

    choice='No';
    while (strcmpi(choice, 'No'))
        h=figure('position',[80 80 1000 600]);
        plot(x,y,'o-');
        hold on
        %ask the user to pick the initial parameters graphically
        title(['Please click on top of the ', num2str(Npeaks), ' peaks (left to right)']);
        [x_click,y_click]=ginput(Npeaks);
        
        for kk=1:Npeaks
                idx_x0=(kk-1)*3+1;
                idx_A=(kk-1)*3+3;
                
                param0(idx_x0)=x_click(kk);
                param0(idx_A)=y_click(kk)*intensity_constant;
        end
        param0_extended=param0;
        param0_extended(length(param0)+1)=existing_spectrum(3);
        param0_extended(length(param0)+2)=existing_spectrum(6);
        
        plot(x, curve_providePDMS(param0_extended,x),'r-','linewidth',2)
        title('Initial fit guess');
        
        choice = questdlg('Is the initial fit good?', ...
            'Good initial fit?','Yes','No','Yes');
        
        close(h);
    end    
end

%%
lb=[1,2,1, 1,2,1,    1,1]; %pos, width, intensity
ub=[max(x),50,max(y)*1e2, max(x),50,max(y)*1e2,    max(y)*1e2,max(y)*1e2];
lb(1:3:6)=param0_extended(1:3:6)-3;
ub(1:3:6)=param0_extended(1:3:6)+3;
% lb(2:3:6)=param0(2:3:6);
% ub(2:3:6)=param0(2:3:6);

options=optimoptions('lsqcurvefit','display','off','algorithm','trust-region-reflective',...
    'MaxFunEvals',1e4, 'maxiter', 1e5, 'TolX', 1e-6, 'TolFun', 1e-6);

[param,fval]=lsqcurvefit(fun,param0_extended,x,y,lb,ub,options);
%%

if plotON
    xfine=linspace(x(1),x(end),200);
    figure(2)
    plot(x,y,'bo-');
    hold on
    % plot(x, curve(param,x),'r-');
    plot(xfine,curve(param,xfine),'r-');
    plot(xfine,curve_providePDMS(param,xfine),'b-');
    plot(param(1:3:6), repmat(max(y),length(param(1:3:6))), 'm*');
    axis([x(1) x(end) 0 max(y)*1.2]);
    hold off
end
peak_positions=param(1:3:6);
end

