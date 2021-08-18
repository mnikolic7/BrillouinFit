function [merged] = merge_unequal(a1,varargin)
%this function merges uneven column arrays into a matrix of size (max leng
% of arrays) x (num arrays). pads with nans
%   author: mnikolic@umd.edu
N=nargin;
L=zeros(1,N);
L(1)=size(a1,1);
for k=2:N
    L(k)=size(varargin{k-1},1);
end

Nrows=max(L);

merged=nan(Nrows,N);
merged(1:L(1),1)=a1;
for k=2:N
    merged(1:L(k),k)=varargin{k-1};
end


end

