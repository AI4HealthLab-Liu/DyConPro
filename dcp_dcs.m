function [dcs,ifc]=dcp_dcs(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.0 release 10/18/2017

ifc=dcp_ifc(x);
[~,~,dcs]=dcp_get_envelope(ifc);
