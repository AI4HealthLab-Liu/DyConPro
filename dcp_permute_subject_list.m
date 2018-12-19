function out_list=dcp_permute_subject_list(x,nsubs)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Inputs:
% 1. x is a list of subject numbers
% 2. nsubs is #subs to include in the sample per experimet
% 
% Output:
% 1. out_list is a matrix of permuted subject ID #'s
% 

col_count=0;
out_list=zeros(nsubs,floor(length(x)/nsubs));
while length(x)>nsubs
    col_count=col_count+1;
    perminds=randperm(length(x),nsubs);
    out_list(:,col_count)=x(perminds);
    x(perminds)=[];
end



end