function [data_matrix] = getDataMatrix(field_name, varargin)
%This function takes in the data structs from function
%calculate_cell_metrics_shapeBF and similar. The data struct is usually
%named "vJan01_2D" or variation of this. This function takes in multiple
%structs and makes a column for each one of them. The first argument is the
%field name that is to be extracted. 
%   author: mnikolic@umd.edu

data_matrix=nan(25,nargin-1);

for n=1:nargin-1
    dStruct=varargin{n};
    M=[dStruct.(field_name)]; %amazing syntax. love it. 
    if length(M)>25 %if M is longer than 25 (extremely rare and unlikely)
        data_matrix2=nan(length(M),nargin-1); %initialize a large enough matrix
        data_matrix2(1:25,:)=data_matrix; %copy the present data into it
        data_matrix=data_matrix2; %and replace the current matrix with the larger one.
    end
    data_matrix(1:length(M),n)=M;
end

end

