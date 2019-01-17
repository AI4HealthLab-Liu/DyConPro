% function x=dcp_coord_table(P,thresh,roivec,roi_name,nodefile)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox

% INPUTS:
% 1. P is the set of parafac factors
% 2. thresh is the zmap threshold
% 3. roivec is a list of rois (see example below)
% 4. roi_name is a string such as 'vmpfc'
% 5. nodefile must be a .txt file extension of a BNV .node type file

% % anterior cingulate labels at 5, 22, 40
% roivec=[];
% roivec=[5 6 22 29 40 46 55 76 79]; %ACCs
% roivec=[51 109]; %vmPFC
roivec=[83 112]; % Insula
roi_name='Insula';
if ~isempty(roivec)==0
    comps=length(P.lambda);
end
if ~isempty(roivec)==1
    comps=1;
end
thresh=2.5;
nodefile='MNINODEFILE.txt';
% 
for fac=[8] %1:comps
    map=P{2}(:,fac)*P{3}(:,fac)';
    Zmap=dcp_zscore_thresh(map);
    if ~isempty(roivec)==1
        ROImap=zeros(size(Zmap));
        ROImap(:,roivec)=Zmap(:,roivec); % This makes a triu of Zmap for ROI coords only
        ROImap=ROImap+ROImap'; % This makes it a symmetric matrix
        save(['comp',int2str(fac),'_',roi_name,'_zmap.edge'],'ROImap','-ascii')
    end
    if ~isempty(roivec)==0
        save(['comp',int2str(fac),'_zmap.edge'],'Zmap','-ascii')
    end
end

% % Generate tables of coordinates and labels
for cnum=[8] %1:comps
    T=readtable(nodefile,'ReadVariableNames',false);
    if ~isempty(roivec)==1
        X=load(['comp',int2str(cnum),'_',roi_name,'_zmap.edge']); 
    end
    if ~isempty(roivec)==0
        X=load(['comp',int2str(cnum),'_zmap.edge']);
    end
    Xbp=X;Xbn=X;Xbp(Xbp<thresh)=0;Xbn(Xbn>-thresh)=0;
    node_strg_pos=sum(Xbp); %node strength is sum of weights
    node_strg_neg=sum(Xbn);
    Xbp(Xbp~=0)=1;Xbn(Xbn~=0)=1;
    node_deg_pos=sum(Xbp); %node degree is sum of suprathreshold edges
    node_deg_neg=sum(Xbn);
    regp=find(node_deg_pos>0);regn=find(node_deg_neg>0);
    NodesPos=[T(regp,:) array2table(node_deg_pos(regp)','VariableNames',{'Deg'}) array2table(node_strg_pos(regp)','VariableNames',{'Strg'})];
    NodesNeg=[T(regn,:) array2table(node_deg_neg(regn)','VariableNames',{'Deg'}) array2table(node_strg_neg(regn)','VariableNames',{'Strg'})];
    NodesPos=sortrows(NodesPos,'Deg','descend');NodesNeg=sortrows(NodesNeg,'Deg','descend');
    NodesPos=table2cell(NodesPos);NodesPos=cell2table(NodesPos,'VariableNames',{'x','y','z','Anat_Label','Deg','Strg'});
    NodesNeg=table2cell(NodesNeg);NodesNeg=cell2table(NodesNeg,'VariableNames',{'x','y','z','Anat_Label','Deg','Strg'});
    if ~isempty(roivec)==1
        writetable(NodesPos,['NodesPos_Comp_',int2str(cnum),'_',roi_name,'.xls']) 
        writetable(NodesNeg,['NodesNeg_Comp_',int2str(cnum),'_',roi_name,'.xls'])
    end
    if ~isempty(roivec)==0
        writetable(NodesPos,['NodesPos_Comp_',int2str(cnum),'.xls']) 
        writetable(NodesNeg,['NodesNeg_Comp_',int2str(cnum),'.xls']) 
    end
    X=triu(X,1);
    [row,col]=find(abs(X)>thresh);
    L=cell(1,length(row));
    for loop1=1:length(row)    
        if ~isempty(roivec)==0
            L{loop1}=[table2cell(T(row(loop1),4)) table2cell(T(row(loop1),1:3)) X(row(loop1),col(loop1)) table2cell(T(col(loop1),1:3)) table2cell(T(col(loop1),4))];  
        end
        if ~isempty(roivec)==1
            L{loop1}=[cell2table(table2cell(T(row(loop1),4)),'VariableNames',{'Anat_A'}) cell2table(table2cell(T(row(loop1),1:3)),'VariableNames',{'x_A','y_A','z_A'}) array2table(X(row(loop1),col(loop1)),'VariableNames',{'Edge_Wt'}) cell2table(table2cell(T(col(loop1),1:3)),'VariableNames',{'x_B','y_B','z_B'}) cell2table(table2cell(T(col(loop1),4)),'VariableNames',{'Anat_B'})];
        end
        idxs(loop1,1)=X(row(loop1),col(loop1));
    end
    for loop2=1:length(row)
        CT(loop2,:)=L{loop2};
    end
    [~,ids]=sort(idxs,'descend');
    sCT=CT(ids,:);
    if ~isempty(roivec)==1
        writetable(sCT,['SuppTable_Comp_',int2str(cnum),'_',roi_name,'_Edges.xls']) 
    end
    if ~isempty(roivec)==0
        writetable(sCT,['SuppTable_Comp_',int2str(cnum),'_Edges.xls']) %,'Delimiter','tab') 
    end
    clear row col sx idxs ids CT sCT X L
end

