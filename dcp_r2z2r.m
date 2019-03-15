function [rzr,rzrt,rzrp]=dcp_r2z2r(Dopts)

% clear,clc

% Dopts.atlasdata
% Dopts.data
% Dopts.dofdata
%     assumes that sub-brik 5 of Dopts.dofdata has dof information
% Dopts.funcdir=[wdir,int2str(sub),'/func'];

% atlas
[~,ATLAS,~,~]=BrikLoad(Dopts.atlasdata);
rois=unique(ATLAS(ATLAS~=0));
[xd,yd,zd]=size(ATLAS);
ATLAS_rs=reshape(ATLAS,1,xd*yd*zd);

% denoised data
[~,V1,~,~]=BrikLoad(Dopts.data);
V1_rs=reshape(V1,xd*yd*zd,size(V1,4))';

% DoF data
if ischar(Dopts.dofdata)
    [~,V2,~,~]=BrikLoad(Dopts.dofdata);
    V2=V2(:,:,:,5);
end
if isnumeric(Dopts.dofdata)
    V2=ones(xd,yd,zd).*Dopts.dofdata;
end

roi_dof=zeros(length(rois),length(rois));
roi_tcs=zeros(length(rois),length(rois));
for roi=1:length(rois)
    roi_dof(:,roi)=size(V1,4)-mean(V2(ATLAS==roi));
    roi_tcs(:,roi)=mean(V1_rs(:,ATLAS_rs==roi),2);
end

[rr,~]=corr(roi_tcs);
rzr=zeros(length(rois),length(rois));
rzrt=zeros(length(rois),length(rois));
rzrp=zeros(length(rois),length(rois));
for loop1=1:length(rois)
    for loop2=1:length(rois)
        dofused=mean([roi_dof(loop1) roi_dof(loop2)]);
        DFnew=size(V1,4)-dofused;
        DForig=size(V1,4);
        zfc=atanh(rr(loop1,loop2));
        zfcT=zfc/(1/(sqrt(DFnew-3)));
        corcorrected=zfcT/(sqrt(DForig-3));
        rzr(loop1,loop2)=tanh(corcorrected);
        rzrt(loop1,loop2)=(rzr(loop1,loop2).*sqrt(DForig-2))./sqrt(1-rzr(loop1,loop2).^2);
        rzrp(loop1,loop2)=2*tcdf(abs(rzrt),(DForig-2),'upper');  
    end
end


% figure();scatter(rr,rnew)
% figure();scatter(rr,rr-rnew)




