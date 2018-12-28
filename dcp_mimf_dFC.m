function [mimf_dfc,ifc_tens,IMF,imf_keep,imf_powspec]=dcp_mimf_dFC(X,w_bot,w_top,f_bot,f_top,Fs)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% DYnamic CONnectivity PROcessing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Inputs:
% 1. X is a time x roi matrix
% 2. w_top & w_bot are bandpass (i.e., full band) for power spectrum calculation
% 3. f_top & f_bot are bandpass for imfs used in mimf_dfc calculations
% 
% Output:
% 1. mimf_dfc is the ifc averaged over imf center_freq bands
% 2. ifc_tens is the ifc for all imf center_freq bands 
%       e.g., (:,:,:,1)=highest freq in centers_keep
% 
% NOTE:
%     This code makes use of the function apitmemd.m to compute IMFs.  To use this function 
%     you must download the code from:
%     

if ~exist('apitmemd.m','file') %|| ~exist('EMDLAB_Plugin','dir')
    error('This function requires the apitmemd.m function from Hemakom et al., 2014');
end

% clear,clc
% % 
% X=randn(167,246);
% % 
% Fs=.5;
% w_bot=.001;
% w_top=.25;
% f_top=.125;
% f_bot=.01;

[td,sigs]=size(X);
[Xmp,mf,~]=dcp_mirror_pad(X);
IMFmp=apitmemd(Xmp,sigs*2,'alpha',1);
IMF=IMFmp(:,:,mf+1:end);
IMF=IMF(:,:,1:td);


[~,nimfs,~]=size(IMF);

imf_power=zeros(sigs,nimfs);
imf_powspec=zeros(sigs,nimfs,td);
for loop1=1:nimfs
    for loop2=1:sigs
    [z,pts]=dcp_mkfreq(squeeze(IMF(loop2,loop1,:)),w_bot,w_top,Fs);
    imf_power(loop2,loop1)=z(pts==max(pts));
    imf_powspec(loop2,loop1,:)=pts;
    end
end
mimf_pow=mean(imf_power);
centers_keep=find(mimf_pow<=f_top & mimf_pow>=f_bot);
imf_keep=IMF(:,centers_keep,:);

ifc_tens=zeros(td,sigs,sigs,length(centers_keep));
for loop1=1:length(centers_keep)
    M=squeeze(imf_keep(:,loop1,:))';
    [ifc]=dcp_ifc(M);
    ifc_tens(:,:,:,loop1)=ifc;
end

mimf_dfc=mean(ifc_tens,4);


