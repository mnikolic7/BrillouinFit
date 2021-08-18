function [data_cell] = getDataCell(field_name, varargin)
%This function takes in the data structs from function
%calculate_cell_metrics_shapeBF and similar. The data struct is usually
%named "vJan01_2D" or variation of this. This function takes in multiple
%structs and makes a cell column for each one of them. The first argument is the
%field name that is to be extracted. 
%   author: mnikolic@umd.edu

data_cell=cell(25,nargin-1);

for n=1:nargin-1
    dStruct=varargin{n};
    
    Nstruct_entries=length(dStruct);
     
    if Nstruct_entries>25 %if M is longer than 25 (extremely rare and unlikely)
        data_cell2=cell(Nstruct_entries,nargin-1); %initialize a large enough matrix
        data_cell2(1:25,:)=data_cell; %copy the present data into it
        data_cell=data_cell2; %and replace the current matrix with the larger one.
    end
    
    for k=1:Nstruct_entries
        field=dStruct(k).(field_name);
        data_cell{k,n}=field;
    end
end

end

