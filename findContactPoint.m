function [bad_curve_idx,contact_points,toe_region_ends,idxCP,idxTRE] = findContactPoint(p_all,f_all)
%This function takes in the interpolated force curves as a matrix of column
%vectors: positions (p_all) and forces (f_all). It calculates the numerical
% derivative of the force. At large distances there is zero slope, and at
% large indentation there is constant negative slope. This function fits an
% function to the first 70% of the derivative curve to find the transition
% between these two plateaus. Then it identifies the contact point as point
% 2*sigma from the center of the. The end of the toe region starts one
% sigma in the other direction from the center of the erf.
% Finally, this function displays the force curve with points identified by
% circles. It will do this for each column of p_all and f_all.
%   mnikolic@umd.edu
N=size(p_all,2);
contact_points=nan(N,1);
toe_region_ends=nan(N,1);
bad_curve_idx=zeros(N,1);
contact_points=zeros(N,1);
toe_region_ends=zeros(N,1);

for k=1:N
    x=p_all(:,k);
    y=f_all(:,k);
    gy=smooth(gradient(y,x(2)-x(1)));
    
    xf=x(1:700);
    gyf=gy(1:700);
    [xData, yData] = prepareCurveData( xf, gyf );
    
    % Set up fittype and options.
    ft = fittype( '(a/2)*(1+erf((x-b)/(c*sqrt(2))))+d;', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [0.004 3e-06 2e-06 -0.008];
    
    % Fit model to data.
    [fitresult, ~] = fit( xData, yData, ft, opts );
    
    plot(x,y);
    hold on
    % plot(fitresult,'r-');
    [~,idxC]=min(abs(x-fitresult.b));
    plot(x(idxC),y(idxC),'ro');
    [~,idxcp]=min(abs(x-(fitresult.b+fitresult.c*2)));
    plot(x(idxcp),y(idxcp),'ro');
    [~,idxtre]=min(abs(x-(fitresult.b-fitresult.c)));
    plot(x(idxtre),y(idxtre),'ro');
    title(['Curve number: ', num2str(k)]);
    hold off
    idxCP(k)=idxcp;
    idxTRE(k)=idxtre;
    answer = questdlg('Is the curve good', ...
        'Curve review', ...
        'Yes','No','No');
    % Handle response
    switch answer
        case 'Yes'
            bad_curve_idx(k) = 0;
        case 'No'
            bad_curve_idx(k) = 1;
            display(['Ok. Curve number: ', num2str(k),' marked as bad']);
    end
    contact_points(k)=x(idxcp);
    toe_region_ends(k)=x(idxtre);
    
end
end

