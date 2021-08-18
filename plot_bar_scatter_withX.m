function plot_bar_scatter_withX(data,x_pos,groups,cols,statsOn,varargin)
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
% all_colors=zeros(8,3);
% % colors={'0b00ec','bd00c0','fb008e','ff0060','ff3d38','ff840d','ffb600'};
% colors={'003f5c','374c80','7a5195','bc5090','ef5675','ff764a','ffa600'};
% for k=1:7
%     all_colors(k,:)=reshape(sscanf((colors{k}).','%2x'),3,[]).'/255;
% end

cmap=colormap('jet');
Ndata=size(data,2);
cmapidx=round(linspace(65,192,Ndata));
all_colors=cmap(cmapidx,:);


kdeON=1;
PlotPointsOn=1;
%%

if nargin==6
    ax=varargin{1};
else
    ax=gca;
end

Ndata=size(data,2);


hold on
for pos_idx=1:Ndata
%         cidx=mod(pos-1,8)+1;
    cidx=pos_idx;
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
    
    if kdeON==1
        col2=col;
        minb=min(data(:,pos_idx),[],'omitnan');
        maxb=max(data(:,pos_idx),[],'omitnan');
        Np=length(data(:,pos_idx));
%         bin_width=(maxb-minb)./Np*12
        bin_width=0.01;
        minb=minb;%-2.5*bin_width;
        maxb=maxb;%+2.5*bin_width;
        bin_centers=6:bin_width:7;
        bin_centers_fine=linspace(minb,maxb,200);
        pdKDE = fitdist(data(:,pos_idx),'kernel','BandWidth',bin_width);
        yKDE=pdf(pdKDE, bin_centers_fine)*bin_width;
        yKDE=0.5*yKDE./max(yKDE(:));
%         text(x_pos(pos_idx),6.13,['n=',num2str(sum(~isnan(data(:,pos_idx))),'%2.d')],'horizontalalignment','center','fontsize',8);
        horz=[-yKDE fliplr(yKDE)];
        vert=[bin_centers_fine, fliplr(bin_centers_fine)];
        hKDE=fill(ax,x_pos(pos_idx)+horz,vert,col2,'edgecolor','none');
        alpha(hKDE,0.75);
        hlkde=plot(ax,x_pos(pos_idx)+[-0.2 0.2], [nanmean(data(:,pos_idx)),nanmean(data(:,pos_idx))],'-','linewidth',1,'color','k');
%         hekde1=plot(ax,x_pos(pos_idx)+[-0.2 0.2], [1 1]*(nanmean(data(:,pos_idx))+nanstd(data(:,pos_idx))),'-','color','k','linewidth',1);
%         hekde2=plot(ax,x_pos(pos_idx)+[-0.2 0.2], [1 1]*(nanmean(data(:,pos_idx))-nanstd(data(:,pos_idx))),'-','color','k','linewidth',1);
        he=errorbar(ax,x_pos(pos_idx), nanmean(data(:,pos_idx)), nanstd(data(:,pos_idx)),'-','color','k');
    end
    
    
     if PlotPointsOn==1
        
        col3=col;
        col3=[1 1 1]*0.25;
        %     col2=0.7*[1 1 1];
        hp=plot(ax,x_pos(pos_idx)-2*0.125+2*0.25*rand(size(data(:,pos_idx))),data(:,pos_idx),'ko','linewidth',0.25,'color',col3,'markerfacecolor',col3,'markersize',1);
%         he=errorbar(ax,x_pos(pos_idx), nanmean(data(:,pos_idx)), nanstd(data(:,pos_idx)),'-','color','k','linewidth',1);
%         hl=plot(ax,x_pos(pos_idx)+[-0.2 0.2], [nanmean(data(:,pos_idx)),nanmean(data(:,pos_idx))],'-','linewidth',3,'color','k');
        %     alpha(hp,transparency);
        %     alpha(he,transparency);
        
        
    end
end
% set(ax,'xtick',x_pos,'xticklabel',groups);
set(ax,'xtick',0:8)
% ylabel('Average Brillouin Shift (GHz)','fontsize',14);
% set(gca,'fontsize',14)
% xlabel(' ','FontSize',15)
% axis([0 Ndata+1 nanmin(data(:))*0.8, nanmax(data(:))*1.2]);
%
% set(ax,'linewidth',1.5);
% set(ax,'fontsize',14);
if statsOn
    [p,tbl,stats] = anova1(data,[],'off');
    %     disp('test = kruskalwallis');
    disp('test = anova1');
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

YMIN=min(data(:));
YMAX=max(data(:));
RANG=YMAX-YMIN;
set(ax,'ylim',[YMIN-0.2*RANG YMAX+0.2*RANG]);
hold off
end

