function plot_bar_scatter_KDE(data,groups,cols,statsOn,varargin)
%PLOT_BAR_SCATTER plots bar graph of data with group names groups
%   author: mnikolic@umd.edu

% color_blue=[0, 0.4470, 0.7410];
% color_red=[0.8500, 0.3250, 0.0980];
% color_orange=[0.9290, 0.6940, 0.1250];
% color_purple=[0.4940, 0.1840, 0.5560];
% color_green=[0.4660, 0.6740, 0.1880];
% color_lightblue=[0.3010, 0.7450, 0.9330];
% color_darkred=[0.6350, 0.0780, 0.1840];
% color_black=[0, 0, 0];
% all_colors=[color_blue;color_red;color_orange;color_purple;color_green;color_lightblue;color_darkred;color_black];
all_colors=zeros(8,3);
% colors={'0b00ec','bd00c0','fb008e','ff0060','ff3d38','ff840d','ffb600'};
colors={'003f5c','374c80','7a5195','bc5090','ef5675','ff764a','ffa600'};
for k=1:7
    all_colors(k,:)=reshape(sscanf((colors{k}).','%2x'),3,[]).'/255;
end
%%

if nargin==5
    ax=varargin{1};
else
    ax=gca;
end

Ndata=size(data,2);


hold on
for pos=1:Ndata
%     cidx=mod(pos-1,8)+1;
        cidx=cols(pos);
    col=all_colors(cidx,:);
%         hb=bar(pos, nanmean(data(:,pos)),'edgecolor','none','facecolor',col,'barwidth',0.7);
%         alpha(hb,0.3);
%     %%%%plot gaussian distribution random around the x position
%         y=data(:,pos);
%         dy=(y-nanmean(y))./nanstd(y);
%         f_y=exp(-abs(dy).^2);
%         xpositions=pos+0.25.*f_y.*(2*rand(size(y))-1);
%         hp=plot(xpositions,data(:,pos),'ko','linewidth',0.25,'color',col,'markerfacecolor',col,'markersize',1);
%     %%%%%%
    col2=col.^(1.5);
%     col2=0.7*[1 1 1];
%     hp=plot(ax,pos-0.125+0.25*rand(size(data(:,pos))),data(:,pos),'ko','linewidth',0.25,'color',col2,'markerfacecolor',col2,'markersize',2);
%     he=errorbar(ax,pos, nanmean(data(:,pos)), nanstd(data(:,pos)),'-','color','k','linewidth',1.5);
%     hl=plot(ax,[pos-0.2 pos+0.2], [nanmean(data(:,pos)),nanmean(data(:,pos))],'-','linewidth',2,'color',col);
    
    minb=min(data(:,pos));
    maxb=max(data(:,pos));
    Np=length(data(:,pos));
    bin_width=(maxb-minb)./Np*10;
    minb=minb-2.5*bin_width;
    maxb=maxb+2.5*bin_width;
    bin_centers=minb:bin_width:maxb;
    pdKDE = fitdist(data(:,pos),'kernel','BandWidth',bin_width);
    yKDE=pdf(pdKDE, bin_centers)*bin_width;
    yKDE=0.25*yKDE./max(yKDE(:));
    
%     text(pos,maxb+2*bin_width,['n=',num2str(sum(~isnan(data(:,pos))),'%2.d')],'horizontalalignment','center');
    text(pos,6.25,['n=',num2str(sum(~isnan(data(:,pos))),'%2.d')],'horizontalalignment','center');
    
%     hKDE=plot(ax,pos-0.3-yKDE*0.2,bin_centers,'-','color',col);
    horz=[-yKDE fliplr(yKDE)];
    vert=[bin_centers, fliplr(bin_centers)];
    hKDE=fill(ax,pos+horz,vert,col,'edgecolor','none');
    alpha(hKDE,0.2);
    
    
    %     alpha(hp,transparency);
    %     alpha(he,transparency);
end
set(ax,'xtick',1:Ndata,'xticklabel',groups);
% ylabel('Average Brillouin Shift (GHz)','fontsize',14);
% set(gca,'fontsize',14)
% xlabel(' ','FontSize',15)
axis([0 Ndata+1 nanmin(data(:))*0.8, nanmax(data(:))*1.2]);
% set(ax,'ylim',[6.1 6.5]);
set(ax,'linewidth',1.5);
set(ax,'fontsize',14);
if statsOn
    [p,tbl,stats] = kruskalwallis(data,[],'off');
    disp('test type = kruskalwallis');
    c = multcompare(stats,'Display','off');
    % comparison_idx=[1  8  11 ];
    comparison_idx=find(c(:,6)<0.05);
    pvals=zeros(1,length(comparison_idx));
    comparison_groups=cell(1,length(comparison_idx));
    
    for K=1:length(comparison_idx)
        pvals(K)=c(comparison_idx(K),6);
        comparison_groups{K}=[c(comparison_idx(K),1),c(comparison_idx(K),2)];
    end
    sigstarOut=sigstar(comparison_groups,pvals,0);
end
hold off
end

