function [radial_image,r_scan,th_scan] = plot_radial_profile_forPCA(x_scan,y_scan,scan_image,bw)
%This function takes in the scan image from a brillouin shift measurement
%and it then returns the image converted to radial (y
%axis) and tangential (horizontal axis). It also returns
%
%   input: x_scan,y_scan (coordinates in um)
%          scan_image: shift image
%          bw: binary image whose boundary will be used to transform the
%          image
%
%   output: radial_image - converted image: y axis is distance from edge to
%   center, x axis is angle (arbitrary start). 
%            r_scan, and th_scan are the coordinates of the transformed
%            image in microns and degrees respectively.
%           r_scan is not exact since r_scan(um) varies across the boundary of
%           the spheroid. So r_scan is just approximate um sized image. 
%           if you want exact r_scan, you should not use this function but
%           just plot the image in spherical coordinates with center at the
%           centroid. This is normalized with respect to the boundary of
%           the spheroid (which is more important).
%
% you can also include teh widht image, but I am thinking shift right now,
% and everywhere where I say shift I really just mean pixel intensity.
%   author: mnikolic@umd.edu
%%
Npts_radial=200;
Npts_tangential=200;
%%
%residual code - identify bw based on scan_image
% bw6p4=scan_image>t;
% bw6p4=imdilate(bw6p4,strel('disk',3));
% bw6p4=imfill(bw6p4,'holes');
% bw6p4=imerode(bw6p4,strel('disk',3));
% bw=bw6p4;

%% enlarge the image (pad with nans) to avoid edge issues
[rows,cols]=size(scan_image);
margin=250;
scan_image2=nan(rows+2*margin,cols+2*margin);
bw2=zeros(rows+2*margin,cols+2*margin);

idx_r=(1:rows)+margin;
idx_c=(1:cols)+margin;
scan_image2(idx_r,idx_c)=scan_image;
bw2(idx_r,idx_c)=bw;

xumpp=x_scan(2)-x_scan(1);
yumpp=y_scan(2)-y_scan(1);

x_scan_pos0=margin*xumpp;
y_scan_pos0=margin*yumpp;
x_scan2=(1:size(scan_image2,2))*xumpp - x_scan_pos0;
y_scan2=(1:size(scan_image2,2))*yumpp - y_scan_pos0;
%% find the centroid
s=regionprops(bw2,'Centroid');
centroid=s.Centroid;
centroid=round(centroid);
%% find the boundaries
B=bwboundaries(bw2);
b=B{1};
xb=b(:,2);
yb=b(:,1);
centroid=round(mean(b));
%%
% imagesc(bw2);
% hold on
% plot(centroid(2),centroid(1),'ro')
% plot(xb,yb,'r-');
% hold off
%%
% imagesc(bw);
%extend each line by 5 um away from the edge y=m*x+b
dy=(centroid(1)-yb);
dx=(centroid(2)-xb);

dy2=-dy./sqrt(dx.^2+dy.^2);
dx2=-dx./sqrt(dx.^2+dy.^2);


delta=0; %um
delta_px=delta/xumpp;
x0=xb+dx2*delta_px;
y0=yb+dy2*delta_px;

% residual code that evenly distributes the theta angles, but it has some 
% boundary case issues with interpolating r from theta values. so those
% need to be fixed before we use it. The current code just interpolates x,y
% values of the boundary position, so that should still work, and is not
% the main issue. 

% r=sqrt((y0-centroid(1)).^2+(x0-centroid(1)).^2);
% theta=atan2(y0-centroid(1),x0-centroid(1));
% theta=unwrap(theta);
% 
% % Npts_radial=length(r);
% % Npts_tangential=length(r);
% theta2=linspace(-pi,pi,Npts_radial); %uniformly distribute angles.
% % theta2=interp1(1:length(theta),theta,linspace(1,length(theta),Npts_radial),'linear');
% r2=interp1(theta,r,theta2,'linear');
% r2=fillmissing(r2,'nearest');
% 
% x0=r2.*cos(theta2)+centroid(2);
% y0=r2.*sin(theta2)+centroid(1);
%%
x0=interp1(x0,linspace(1,length(x0),Npts_radial),'nearest');
y0=interp1(y0,linspace(1,length(y0),Npts_radial),'nearest');

x0=round(x0);
y0=round(y0);

%
%
c=nan(Npts_radial,Npts_tangential);
% hold on
for k=1:length(x0)
    x_line_px=[x0(k) centroid(2)];
    y_line_px=[y0(k) centroid(1)];
%     x_line=x_scan2(x_line_px);
%     y_line=y_scan2(y_line_px);
    
%     c(:,k)=improfile(x_scan2,y_scan2,scan_image2,x_line,y_line,Npts_radial,'bicubic');
    c(:,k)=improfile(scan_image2,x_line_px,y_line_px,Npts_radial,'bicubic');
%     plot((1:100)*0.5,c(:,k),'-');
%     plot(x_line,y_line,'-');
%     pause(0.1)
end
% hold off
%%
radial_image=c;
th_scan=linspace(0,360,Npts_tangential);
avg_r_distance=mean(sqrt((x0-centroid(2)).^2+(y0-centroid(1)).^2));
r_scan=xumpp*avg_r_distance*linspace(0,1,Npts_radial);
% r_scan=r_scan; %the spheroid boundary starts at r=0.

%%
% imagesc(th_scan,r_scan,radial_image);
% 
end

