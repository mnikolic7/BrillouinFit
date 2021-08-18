function [RP,RPerr,RP_ncells] = average_radial_profiles(data_cell)
%This function takes in the data cell from getDataCell.m function. This
%should have n columns (one for each experiment), and a variable number of
%rows (one for each spheroid/cell imaged). This function reads the contents
%of each cell (which should be 200by1 vectors of radial distance) and
%averages them together for all cells in the given experiment (column of
%the data_cell). It omits the empty cells {0x0}. In the end it returns the
%average profile+error across cells.
%   author:mnikolic@umd.edu

Nexp=size(data_cell,2);
RP=nan(200,Nexp);
RPerr=nan(200,Nexp);
RP_ncells=nan(1,Nexp);
for jj=1:4
    curr_rp_all=nan(200,1);
    ii=1;
    for k=1:25
        curr_rp=data_cell{k,jj};
        if isempty(curr_rp)
            continue;
        end
        curr_rp_all(:,ii)=curr_rp;
        ii=ii+1;
    end
    RP(:,jj)=nanmean(curr_rp_all,2);
    n=length(curr_rp_all(1,:));
    RPerr(:,jj)=nanstd(curr_rp_all,[],2)./sqrt(n);
    RP_ncells(1,jj)=n;
end
    
end

