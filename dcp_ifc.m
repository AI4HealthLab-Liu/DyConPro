function [ifc]=dcp_ifc(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018

[row,col]=size(x);
[x1,xfr,~]=dcp_mirror_pad(x);

phx=angle(hilbert(x1));
phx=phx(xfr+1:xfr+row,:);
ifc=zeros(row,col,col);

for slice=1:col
    ifc(:,:,slice)=bsxfun(@minus,phx(:,slice),phx);
end
ifc=cos(ifc);

end
