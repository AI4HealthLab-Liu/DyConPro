function [glob_sync,kuramoto]=dcp_kuramoto(X)

% X is a 3D dFC tensor from a single subject
% kuramoto is global mean dFC time series

[tt,~,~]=size(X);
kuramoto=zeros(tt,1);
glob_sync=zeros(tt,1);

for loop1=1:tt
    temp=triu(squeeze(X(loop1,:,:)),1);
    temp=temp(find(~tril(ones(size(temp)))));
    glob_sync(loop1)=mean(temp);
    kuramoto(loop1)=std(temp);
end

end