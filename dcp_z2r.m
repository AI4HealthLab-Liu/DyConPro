function z2r=dcp_z2r(X)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 

[td,rois]=size(X);
z2r=zeros(td,roi);
for loop=1:rois
    z2r(:,loop)=(exp(2*X(:,loop))-1)./(exp(2*X(:,loop))+1);
end

end