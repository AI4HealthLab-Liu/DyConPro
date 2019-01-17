% function [booboo]=dcp_parafac_collage(Factors,nodefile,edgefile,outfile,thresh)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox

% Inputs:
% 1. Factors is a cell array of factors from CMTF with orthonormal grouping basis
% 2. nodefile is a BNV node file for example the MNINODEFILE for the Craddock200 parcellation
% 3. edgefile is a BNV edge file corresponding to the 'matrix' representation of your brain network
%       !!! This is no longer a required input. Edge file is computed
%           internally from Factors matrices
% 4. outfile name is the name of all outputted files and figures

% Example inputs:
%     nodefile='MNINODEFILE.node';
%     brainmeshfile='BrainMesh_Ch2withCerebellum.nv';
%     outfile='COMO4_R25';

% Outputs:
% There aren't really any, but it makes and saves figures for you
% This script/function is kind of slow, but it goes as fast as BrainNet Viewer will let you go
% 

Factors=P15;
comps=length(Factors.lambda);
nodefile='MNINODEFILE.node';
brainmeshfile='BrainMesh_Ch2withCerebellum.nv';
outfile='COMO1_mo24_R15';
thresh=.0695;thline=thresh.*ones(1,225);

% % make spatial factor loading matrices, figs & edge files
for loop2=1:comps
    X=Factors{2}(:,loop2)*Factors{3}(:,loop2)'; % THIS IS FOR PARAFAC
    XZ=dcp_zscore_thresh(X); % matrix is not thresh'd for image
    XZM(:,:,loop2)=XZ;
    THIS=figure();imagesc(XZ);colormap jet;colorbar;caxis([-3 3])
    saveas(THIS,[outfile,'_comp_',int2str(loop2),'.jpg'])
    XZ(abs(XZ)<2.5)=0; % Matrix is thresh'd for edge file
    save([outfile,'_comp_',int2str(loop2),'.edge'],'XZ','-ascii')
end
close all

% make brain map figures .jpg
for facs=1:comps
    BrainNet_MapCfg(brainmeshfile,nodefile,[outfile,'_comp_',int2str(facs),'.edge'],'BNV_Config_LeftLat.mat',[outfile,'_LeftLat_comp_',int2str(facs),'.jpg']);
    BrainNet_MapCfg(brainmeshfile,nodefile,[outfile,'_comp_',int2str(facs),'.edge'],'BNV_Config_RightLat.mat',[outfile,'_RightLat_comp_',int2str(facs),'.jpg']);
    BrainNet_MapCfg(brainmeshfile,nodefile,[outfile,'_comp_',int2str(facs),'.edge'],'BNV_Config_TopAxial.mat',[outfile,'_TopAxial_comp_',int2str(facs),'.jpg']);
    BrainNet_MapCfg(brainmeshfile,nodefile,[outfile,'_comp_',int2str(facs),'.edge'],'BNV_Config_BottomAxial.mat',[outfile,'_BottomAxial_comp_',int2str(facs),'.jpg']);
    BrainNet_MapCfg(brainmeshfile,nodefile,[outfile,'_comp_',int2str(facs),'.edge'],'BNV_Config_PostCoronal.mat',[outfile,'_PostCoronal_comp_',int2str(facs),'.jpg']);
    BrainNet_MapCfg(brainmeshfile,nodefile,[outfile,'_comp_',int2str(facs),'.edge'],'BNV_Config_AntCoronal.mat',[outfile,'_AntCoronal_comp_',int2str(facs),'.jpg']);
end
close(BrainNet)


for facs=1:comps
% Component Brain Maps
    THIS1=figure();
    subplot(3,3,1)
    imshow([outfile,'_LeftLat_comp_',int2str(facs),'.jpg']);
    title('Left Lateral')
    subplot(3,3,2)
    imshow([outfile,'_TopAxial_comp_',int2str(facs),'.jpg']);
    title('Top Axial')
    subplot(3,3,3)
    imshow([outfile,'_RightLat_comp_',int2str(facs),'.jpg']);
    title('Right Lateral')
    subplot(3,3,4)
    imshow([outfile,'_AntCoronal_comp_',int2str(facs),'.jpg']);
    title('Anterior Coronal')
    subplot(3,3,5)
    imshow([outfile,'_BottomAxial_comp_',int2str(facs),'.jpg']);
    title('Bottom Axial')
    subplot(3,3,6)
    imshow([outfile,'_PostCoronal_comp_',int2str(facs),'.jpg']);
    title('Posterior Coronal')

% Component Spatial Matrices Plots  
    subplot(3,3,7)
%     imshow([outfile,'_comp_',int2str(facs),'.jpg'])
    xb=imagesc(XZM(:,:,facs));
    colormap jet;hb=colorbar('Position',[.3395 .1082 .0119 .2],'FontSize',6);caxis([-3 3])
    xlabel('Brain Region');ylabel('Brain Region')
    title('Spatial Loadings')
    set(gca,'FontSize',8)
 
% Component Participant Loadings Plots
    subplot(3,3,8)
    bar(Factors{4}(:,facs));xlim([0 17]);box off
    title('Participant Loadings')
    xlabel('Participant #');ylabel('Factor Loading')
    set(gca,'FontSize',8)

% Component Time Signature Plots
    subplot(3,3,9)
    plot(thline,'k--');hold on
    plot(Factors{1}(:,facs));xlim([0 225]);ylim([min(min(Factors{1})) max(max(Factors{1}))]);box off
    title(['Comp ',int2str(facs),' Time Signature']);xlabel('Time (TRs)');ylabel('Synchronization Strength')
    set(gca,'FontSize',8)

% Figure settings and save THIS1
    set(THIS1,'position',[279 130 1077 852])   
    set(THIS1,'Color',[1 1 1])
    suptitle(['Component ',int2str(facs),': Network & Factor Loadings'])
    set(xb.Parent,'Position',[.1200 .1082 .2134 .2055])
    saveas(THIS1,[outfile,'_BrainCollage_comp_',int2str(facs),'.jpg'])
end

