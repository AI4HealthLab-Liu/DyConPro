function varargout = fc_snapshot(varargin)
% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% FC_SNAPSHOT MATLAB code for fc_snapshot.fig
% 
%     Usage: fc_snapshot(input.struct)
%     
%     Input is a structure with the following fields:
%         1. input.dfc is a dfc tensor or multilayer graph; should be symmetric
%         2. input.coords is a nodes x 3 matrix of node coordinates
%         3. input.nodelabels is a cell array (actually a matrix of strings) of node labels
%         4. input.icnid is a numerical vector indicating network membership
%         5. input.icnlabels is a cell array of network labels for networks in icnid
% 
%      FC_SNAPSHOT, by itself, creates a new FC_SNAPSHOT or raises the existing
%      singleton*.
%
%      H = FC_SNAPSHOT returns the handle to a new FC_SNAPSHOT or the handle to
%      the existing singleton*.
%
%      FC_SNAPSHOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FC_SNAPSHOT.M with the given input arguments.
%
%      FC_SNAPSHOT('Property','Value',...) creates a new FC_SNAPSHOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fc_snapshot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fc_snapshot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fc_snapshot

% Last Modified by GUIDE v2.5 08-May-2018 23:16:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fc_snapshot_OpeningFcn, ...
                   'gui_OutputFcn',  @fc_snapshot_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
% if ~isempty(varargin)
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end
    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
% % % else
% % %     error('Error 1: You need to input some data.')
% % % end
% End initialization code - DO NOT EDIT


% --- Executes just before fc_snapshot is made visible.
function fc_snapshot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fc_snapshot (see VARARGIN)
global circ_thresh

wsvars=evalin('base','who');
if isempty(wsvars)
    wsvars='none';
end
handles.popupmenu2.String=wsvars;

handles.time_editbox.String='1';
handles.thresh_editbox.String='.5';

set(handles.thresh_editbox,'String',.5);
circ_thresh=str2num(get(handles.thresh_editbox,'String')); %#ok<*ST2NM>

% Choose default command line output for fc_snapshot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


%{
% Initialize the input.struct data as global dcp_data for viewing.
global dcp_data dcp_data_orig tframe length_data dimsN circ_thresh G V X Y Z sagslice axslice corslice Vflip
global dcp_coords dcp_nodelabels_cell dcp_icnid dcp_icnlabels dcp_nodelabels dcp_coords_cell dcp_coords_vox
if ~isempty(varargin)
    dcp_data=varargin{1}.dfc;
    dcp_data_orig=dcp_data;
    length_data=size(dcp_data,1);
    dimsN=ndims(dcp_data);
    if dimsN~=3 || length_data<2
        error('Error 2: You need to input a dynamic FC tensor.')
    end
    
    if isfield(varargin{1},'nodelabels')
        dcp_nodelabels=varargin{1}.nodelabels;
        if ~iscell(varargin{1}.nodelabels)
            for loop1=1:size(varargin{1}.coords,1)
%                 dcp_nodelabels_cell{loop1}=[int2str(loop1),'. ',dcp_nodelabels(loop1,:)];
                dcp_nodelabels_cell{loop1}=dcp_nodelabels(loop1,:);
            end
        else
            for loop1=1:size(varargin{1}.coords,1)
%                 dcp_nodelabels_cell{loop1}=[int2str(loop1),'. ',varargin{1}.nodelabels{loop1}];
                dcp_nodelabels_cell{loop1}=varargin{1}.nodelabels{loop1};
            end
        end
        set(handles.uitable1,'Data',dcp_nodelabels_cell');
    end
    
    if isfield(varargin{1},'coords') && ~isempty(varargin{1}.coords)
        dcp_coords=varargin{1}.coords;
        dcp_coords_x=mat2cell(dcp_coords(:,1),ones(1,size(dcp_coords,1)));
        dcp_coords_y=mat2cell(dcp_coords(:,2),ones(1,size(dcp_coords,1)));
        dcp_coords_z=mat2cell(dcp_coords(:,3),ones(1,size(dcp_coords,1)));
        set(handles.uitable1,'Data',[dcp_nodelabels_cell' dcp_coords_x dcp_coords_y dcp_coords_z]);
    end
    
    if isfield(varargin{1},'icnid')
        dcp_icnid=varargin{1}.icnid;
    end
    
    if isfield(varargin{1},'icnlabels')
        dcp_icnlabels=varargin{1}.icnlabels;
    end

%     keyboard
%%%%%%%%%%%%% Load MNI T1 for Glass Brain Contours %%%%%%%%%%%%%%%%%%%%
    [~,V,~,~]=BrikLoad('MNI152_T1_2mm+tlrc');
    [Y,X,Z]=meshgrid(1:109,1:91,1:91);
    
    Vflip=V; %flip(V,2);
    
    sagslice=45;axslice=45;corslice=54;
    
    dcp_coords_vox(:,2)=dcp_coords(:,2)+127; % Y
    dcp_coords_vox(:,3)=dcp_coords(:,3)+79; % Z
    dcp_coords_vox(:,1)=dcp_coords(:,1)+90; % X

    % Generate Sagittal slice images
    contourslice(handles.axes3,Y,X,Z,Vflip,[],45,[],3)
    set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
    daspect(handles.axes3,[1 1 1]);
    
    % Generate Axial slice images
    contourslice(handles.axes4,Y,X,Z,Vflip,[],[],45,3)
    set(handles.axes4,'view',[-90 90],'XLim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
    daspect(handles.axes4,[1 1 1]);
    
    % Generate Coronal slice images
    contourslice(handles.axes5,Y,X,Z,Vflip,54,[],[],3)
    set(handles.axes5,'view',[-90 0],'XLim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
    daspect(handles.axes5,[1 1 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
tframe=1;
set(handles.thresh_editbox,'String',.5);
circ_thresh=str2num(get(handles.thresh_editbox,'String')); %#ok<*ST2NM>

% This sets up the initial plot - only do when we are invisible
% so window can get raised using fc_snapshot.
if strcmp(get(hObject,'Visible'),'off') && dimsN==3
    set(handles.time_editbox,'String',1)
    axes(handles.axes1);cla
    imagesc(squeeze(dcp_data(1,:,:)));colorbar;colormap('jet')
%     set(gca,'YTick',0:length_data-1);
%     yticklabels(dcp_nodelabels);
%     set(gca,'FontSize',3);ytickangle(35);
    dcp_data_graph=squeeze(dcp_data(1,:,:));
    dcp_data_graph(dcp_data_graph<circ_thresh)=0;
    G=graph(dcp_data_graph,'upper');
    axes(handles.axes2);cla
    plot(G,'Layout','circle')
end
%}

% --- Outputs from this function are returned to the command line.
function varargout = fc_snapshot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;


% --- Executes on button press in tplusbutton.
function tplusbutton_Callback(hObject, eventdata, handles)
% hObject    handle to tplusbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_data tframe length_data circ_thresh G dcp_coords_vox X Y Z Vflip dcp_data_orig nodewt
global sagslice axslice corslice isdfc dfc_td dcp_nodelabels_cell dcp_coords_x dcp_coords_y dcp_coords_z
if ~isempty(dcp_data) && isdfc==1
    length_data=dfc_td;
    axes(handles.axes1);cla
    if tframe==length_data
        tframe=1;
    else
        tframe=tframe+1;
    end
    set(handles.time_editbox,'String',tframe)
    circ_thresh=str2num(get(handles.thresh_editbox,'String'));
    imagesc(squeeze(dcp_data(tframe,:,:)));colorbar
    % set(gca,'YTick',0:length_data-1);
    % yticklabels(dcp_nodelabels);
    % set(gca,'FontSize',3);ytickangle(35);
    dcp_data_graph=squeeze(dcp_data(tframe,:,:));
    dcp_data_graph(dcp_data_graph<circ_thresh)=0;
    G=graph(dcp_data_graph,'upper');
    axes(handles.axes2);cla
    plot(G,'Layout','circle')

    %Update sagittal connmap t+
%     dcp_data_graph(dcp_data_graph==tril(dcp_data_graph))=0;
    [I,J]=find(dcp_data_graph>=circ_thresh);
    if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        cla(handles.axes3)
        contourslice(handles.axes3,Y,X,Z,Vflip,[],sagslice,[],3);
        set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
        hold(handles.axes3,'on')
        for loop1=1:size(FF,1)
            plot3(handles.axes3,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
    end

    if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        cla(handles.axes4)
        contourslice(handles.axes4,Y,X,Z,Vflip,[],[],axslice,3);
        set(handles.axes4,'view',[-90 90],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
        hold(handles.axes4,'on')
        for loop1=1:size(FF,1)
            plot3(handles.axes4,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
    end
    
    if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        cla(handles.axes5)
        contourslice(handles.axes5,Y,X,Z,Vflip,corslice,[],[],3);
        set(handles.axes5,'view',[-90 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
        hold(handles.axes5,'on')
        for loop1=1:size(FF,1)
            plot3(handles.axes5,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
    end
    
%     if ~isempty(dcp_nodelabels_cell) && ~isempty(dcp_data)
%         dcp_data_thr=dcp_data;
%         dcp_data_thr(dcp_data_thr<circ_thresh)=0;
%         nodewt=cellstr(num2str(sum(squeeze(dcp_data_thr(tframe,:,:)),2)));
%         set(handles.uitable1,'Data',[dcp_nodelabels_cell dcp_coords_x dcp_coords_y dcp_coords_z cell(size(dcp_coords_x)) nodewt]);
%     end
    
end

% --- Executes on button press in tminusbutton.
function tminusbutton_Callback(hObject, eventdata, handles)
% hObject    handle to tminusbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_data tframe length_data circ_thresh G X Y Z dcp_coords_vox Vflip dcp_data_orig nodewt
global sagslice axslice corslice isdfc dfc_td dcp_nodelabels_cell dcp_coords_x dcp_coords_y dcp_coords_z
if ~isempty(dcp_data) && isdfc==1
    length_data=dfc_td;
    axes(handles.axes1);cla
    if tframe==1 || tframe==0
        tframe=length_data;
    else
        tframe=tframe-1;
    end
    set(handles.time_editbox,'String',tframe)
    circ_thresh=str2num(get(handles.thresh_editbox,'String'));
    imagesc(squeeze(dcp_data(tframe,:,:)));colorbar
    % set(gca,'YTick',0:length_data-1);
    % yticklabels(dcp_nodelabels);
    % set(gca,'FontSize',3);ytickangle(35);
    dcp_data_graph=squeeze(dcp_data(tframe,:,:));
    dcp_data_graph(dcp_data_graph<circ_thresh)=0;
    G=graph(dcp_data_graph,'upper');
    axes(handles.axes2);cla
    plot(G,'Layout','circle')

    %Update sagittal connmap t-
%     dcp_data_graph(dcp_data_graph==tril(dcp_data_graph))=0;
    [I,J]=find(dcp_data_graph>=circ_thresh);
    if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        cla(handles.axes3)
        contourslice(handles.axes3,Y,X,Z,Vflip,[],sagslice,[],3);
        set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
        hold(handles.axes3,'on')
        for loop1=1:size(FF,1)
            plot3(handles.axes3,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
    end

    if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        cla(handles.axes4)
        contourslice(handles.axes4,Y,X,Z,Vflip,[],[],axslice,3);
        set(handles.axes4,'view',[-90 90],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
        hold(handles.axes4,'on')
        for loop1=1:size(FF,1)
            plot3(handles.axes4,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
    end
    
    if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        cla(handles.axes5)
        contourslice(handles.axes5,Y,X,Z,Vflip,corslice,[],[],3);
        set(handles.axes5,'view',[-90 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
        hold(handles.axes5,'on')
        for loop1=1:size(FF,1)
            plot3(handles.axes5,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
    end
    
%     if ~isempty(dcp_nodelabels_cell) && ~isempty(dcp_data)
%         dcp_data_thr=dcp_data;
%         dcp_data_thr(dcp_data_thr<circ_thresh)=0;
%         nodewt=cellstr(num2str(sum(squeeze(dcp_data_thr(tframe,:,:)),2)));
%         set(handles.uitable1,'Data',[dcp_nodelabels_cell dcp_coords_x dcp_coords_y dcp_coords_z cell(size(dcp_coords_x)) nodewt]);
%         dcp_data=dcp_data_orig; 
%     end
    
end

% --- Executes on button press in gotobutton.
function gotobutton_Callback(hObject, eventdata, handles)
% hObject    handle to gotobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global tframe dcp_data length_data circ_thresh G sagslice X Y Z dcp_coords_vox Vflip nodewt dcp_data_orig
global axslice corslice isdfc dfc_td dcp_nodelabels_cell dcp_coords_x dcp_coords_y dcp_coords_z

if ~isempty(dcp_data) && isdfc==1
    if ~isempty(str2num(get(handles.time_editbox,'String')))
        length_data=dfc_td;
        if str2num(get(handles.time_editbox,'String'))<=length_data && str2num(get(handles.time_editbox,'String'))>0 
            axes(handles.axes1);cla
            if tframe~=str2num(get(handles.time_editbox,'String'))
                tframe=str2num(get(handles.time_editbox,'String'));
            end
            imagesc(squeeze(dcp_data(tframe,:,:)));colorbar
    %         set(gca,'YTick',0:length_data-1);
    %         yticklabels(dcp_nodelabels);
    %         set(gca,'FontSize',3);ytickangle(35);
            circ_thresh=str2num(get(handles.thresh_editbox,'String'));
            dcp_data_graph=squeeze(dcp_data(tframe,:,:));
            dcp_data_graph(dcp_data_graph<circ_thresh)=0;
            G=graph(dcp_data_graph,'upper');
            axes(handles.axes2);cla
            plot(G,'Layout','circle')

        elseif str2num(get(handles.time_editbox,'String'))>length_data
            set(handles.time_editbox,'String',tframe)  
        elseif str2num(get(handles.time_editbox,'String'))<0
            set(handles.time_editbox,'String',tframe)  
        end

        dcp_data_graph=squeeze(dcp_data(tframe,:,:));
%         dcp_data_graph(dcp_data_graph==tril(dcp_data_graph))=0;
        [I,J]=find(dcp_data_graph>=circ_thresh);
        if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
            FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
            cla(handles.axes3)
            contourslice(handles.axes3,Y,X,Z,Vflip,[],sagslice,[],3);
            set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
            hold(handles.axes3,'on')
            for loop1=1:size(FF,1)
                plot3(handles.axes3,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
            end
    %     daspect(handles.axes3,[1 1 1]);
        end
        
        if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
            [I,J]=find(dcp_data_graph>=circ_thresh);
            FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
            cla(handles.axes4)
            contourslice(handles.axes4,Y,X,Z,Vflip,[],[],axslice,3);
            set(handles.axes4,'view',[-90 90],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
            hold(handles.axes4,'on')
            for loop1=1:size(FF,1)
                plot3(handles.axes4,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
            end
        %     daspect(handles.axes4,[1 1 1]);
        end

        if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
            [I,J]=find(dcp_data_graph>=circ_thresh);
            FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
            cla(handles.axes5)
            contourslice(handles.axes5,Y,X,Z,Vflip,corslice,[],[],3);
            set(handles.axes5,'view',[-90 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
            hold(handles.axes5,'on')
            for loop1=1:size(FF,1)
                plot3(handles.axes5,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
            end
        %     daspect(handles.axes5,[1 1 1]);
        end
% % % % % % % % % % %         NOT SURE
%         if ~isempty(dcp_data) && ~isempty(dcp_nodelabels_cell)
%             dcp_data_thr=dcp_data;
%             dcp_data_thr(dcp_data_thr<circ_thresh)=0;
%             nodewt=cellstr(num2str(sum(squeeze(dcp_data_thr(tframe,:,:)),2)));
%             set(handles.uitable1,'Data',[dcp_nodelabels_cell dcp_coords_x dcp_coords_y dcp_coords_z cell(size(dcp_coords_x)) nodewt]);
%         else
%             set(handles.time_editbox,'String',tframe)  
%         end
% % % % % % % % % % % 

    end
elseif ~isempty(dcp_data) && isdfc==0
    circ_thresh=str2num(get(handles.thresh_editbox,'String'));
    dcp_data_graph=dcp_data;
    dcp_data_graph(dcp_data_graph<circ_thresh)=0;
    G=graph(dcp_data_graph,'upper');
    axes(handles.axes2);cla
    plot(handles.axes2,G,'Layout','circle')
    
% % % % % % % % % % % % % % % %     
%     dcp_data_graph(dcp_data_graph==tril(dcp_data_graph))=0;
    [I,J]=find(dcp_data_graph>=circ_thresh);
    if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        cla(handles.axes3)
        contourslice(handles.axes3,Y,X,Z,Vflip,[],sagslice,[],3);
        set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
        hold(handles.axes3,'on')
        for loop1=1:size(FF,1)
            plot3(handles.axes3,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
%     daspect(handles.axes3,[1 1 1]);
    end

    if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        cla(handles.axes4)
        contourslice(handles.axes4,Y,X,Z,Vflip,[],[],axslice,3);
        set(handles.axes4,'view',[-90 90],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
        hold(handles.axes4,'on')
        for loop1=1:size(FF,1)
            plot3(handles.axes4,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
    %     daspect(handles.axes4,[1 1 1]);
    end

    if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        cla(handles.axes5)
        contourslice(handles.axes5,Y,X,Z,Vflip,corslice,[],[],3);
        set(handles.axes5,'view',[-90 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
        hold(handles.axes5,'on')
        for loop1=1:size(FF,1)
            plot3(handles.axes5,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
    %     daspect(handles.axes5,[1 1 1]);
    end
% % % % % % % % % % % % % % % % % % % % % % %         
end

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end


% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.fc_snapshot_main)


% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.fc_snapshot_main,'Name') '?'],...
                     ['Close ' get(handles.fc_snapshot_main,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end
delete(handles.fc_snapshot_main)


function time_editbox_Callback(hObject, eventdata, handles)
% hObject    handle to time_editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of time_editbox as text
%        str2double(get(hObject,'String')) returns contents of time_editbox as a double


% --- Executes during object creation, after setting all properties.
function time_editbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to time_editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function thresh_editbox_Callback(hObject, eventdata, handles)
% hObject    handle to thresh_editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh_editbox as text
%        str2double(get(hObject,'String')) returns contents of thresh_editbox as a double


% --- Executes during object creation, after setting all properties.
function thresh_editbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh_editbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% function SaveXYFigBut_Callback(hObject, eventdata, handles)
% % hObject    handle to SaveXYFigBut (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % % % % % AMAZING CODE BELOW !!!!!!
% nfg=figure;this=copyobj(handles.axes1,nfg);set(this,'position','default')
% % nfg=figure;this=copyobj(handles.axes1,nfg);set(this,'position','default');set(nfg,'position',[136 561 1500 366]);


% --- Executes when fc_snapshot_main is resized.
function fc_snapshot_main_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to fc_snapshot_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close fc_snapshot_main.
function fc_snapshot_main_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to fc_snapshot_main (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
clear -global
delete(hObject);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_data tframe
if ndims(dcp_data)>2
    nfg1=figure;imagesc(squeeze(dcp_data(tframe,:,:)));colorbar;colormap('jet')
elseif ndims(dcp_data)==2
    nfg1=figure;imagesc(dcp_data(:,:));colorbar;colormap('jet')    
end
% set(gca,'YTick',0:length_data-1);
% yticklabels(dcp_nodelabels);
% set(gca,'FontSize',3);ytickangle(35);
% keyboard

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_data circ_thresh tframe G
circ_thresh=str2num(get(handles.thresh_editbox,'String'));
dcp_data_graph=squeeze(dcp_data(tframe,:,:));
dcp_data_graph(dcp_data_graph<circ_thresh)=0;
G=graph(dcp_data_graph,'upper');
nfg2=figure;
plot(G,'Layout','circle')
% keyboard


% --------------------------------------------------------------------
function main_menu_menu_Callback(hObject, eventdata, handles)
% hObject    handle to main_menu_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global G X Y Z V sagslice axslice corslice dcp_coords_vox dcp_data tframe circ_thresh Vflip

if ~isempty(Vflip)
    dcp_data_graph=squeeze(dcp_data(tframe,:,:));
    dcp_data_graph(dcp_data_graph<circ_thresh)=0;
    dcp_data_graph(dcp_data_graph==tril(dcp_data_graph))=0;
    figure();
    subplot(1,3,1);%contourslice(Y,X,Z,Vflip,[],sagslice,[]);xlim([1 109]);zlim([1 78]);view([0 0]);
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        contourslice(Y,X,Z,Vflip,[],sagslice,[],3);
        view([0 0]);xlim([1 109]);ylim([1 91]);zlim([1 81]);
        set(gca,'XTick',[],'YTick',[],'ZTick',[]);
        hold on
        for loop1=1:size(FF,1)
            plot3([FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
%         daspect(gca,[1 1 1]);

    subplot(1,3,2);%contourslice(Y,X,Z,Vflip,[],[],axslice);ylim([1 91]);xlim([1 109]);view([-90 90]);
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        contourslice(Y,X,Z,Vflip,[],[],axslice,3);
        view([-90 90]);xlim([1 109]);ylim([1 91]);zlim([1 81]);
        set(gca,'XTick',[],'YTick',[],'ZTick',[]);
        hold on
        for loop1=1:size(FF,1)
            plot3([FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
%         daspect(gca,[1 1 1]);

    subplot(1,3,3);%contourslice(Y,X,Z,Vflip,corslice,[],[]);ylim([1 91]);zlim([1 81]);view([-90 0]);
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        contourslice(Y,X,Z,Vflip,corslice,[],[],3);
        view([-90 0]);xlim([1 109]);ylim([1 91]);zlim([1 81]);
        set(gca,'XTick',[],'YTick',[],'ZTick',[]);
        hold on
        for loop1=1:size(FF,1)
            plot3([FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
%         daspect(gca,[1 1 1]);
end
    

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X Y Z V Vflip dcp_data axslice corslice sagslice dcp_coords_vox circ_thresh tframe
if ~isempty(Vflip)
    figure()
    for loop1=1:2:size(V,3)
        contourslice(Y,X,Z,Vflip,loop1,loop1,loop1,2)
    end
    ylim([1 91]);
    xlim([1 109]);
    zlim([1 91]);
    % view([-90 90]);
    daspect(gca,[1 1 1]);
    % cmap=colormap('white');
    % cscla=1:-1/31:0;cscl=[cscla fliplr(cscla)];
    % cmap=bsxfun(@times,cmap,cscl');
    % colormap(cmap);
    colormap(gray(100));
    view(3)
    if ~isempty(dcp_data) && ~isempty(dcp_coords_vox)
        dcp_data_graph=squeeze(dcp_data(tframe,:,:));
        dcp_data_graph(dcp_data_graph<circ_thresh)=0;
        dcp_data_graph(dcp_data_graph==tril(dcp_data_graph))=0;
        [I,J]=find(dcp_data_graph>=circ_thresh);
        FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
        % contourslice(Y,X,Z,Vflip,[],sagslice,[],3);
        % view([0 0]);xlim([1 109]);ylim([1 91]);zlim([1 81]);
        % set(gca,'XTick',[],'YTick',[],'ZTick',[]);
        hold on
        for loop1=1:size(FF,1)
            plot3([FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
        end
    end
end


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over tplusbutton.
function tplusbutton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to tplusbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function import_data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to import_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function data_processing_menu_Callback(hObject, eventdata, handles)
% hObject    handle to data_processing_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function reset_gui_menu_Callback(hObject, eventdata, handles)
% hObject    handle to reset_gui_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fc_snapshot_resetgui


% --------------------------------------------------------------------
function preprocess_menu_Callback(hObject, eventdata, handles)
% hObject    handle to preprocess_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to fc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function plugin_menu_Callback(hObject, eventdata, handles)
% hObject    handle to plugin_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function parafac_menu_Callback(hObject, eventdata, handles)
% hObject    handle to parafac_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cmtf_menu_Callback(hObject, eventdata, handles)
% hObject    handle to cmtf_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function graph_theory_menu_Callback(hObject, eventdata, handles)
% hObject    handle to graph_theory_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function postprocess_menu_Callback(hObject, eventdata, handles)
% hObject    handle to postprocess_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function decom_dfc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to decom_dfc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function edge_stats_menu_Callback(hObject, eventdata, handles)
% hObject    handle to edge_stats_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function net_stats_menu_Callback(hObject, eventdata, handles)
% hObject    handle to net_stats_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ncnc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to ncnc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function decomp_pca_menu_Callback(hObject, eventdata, handles)
% hObject    handle to decomp_pca_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function decomp_svd_meu_Callback(hObject, eventdata, handles)
% hObject    handle to decomp_svd_meu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function decomp_kmeans_menu_Callback(hObject, eventdata, handles)
% hObject    handle to decomp_kmeans_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function sig_thresh_menu_Callback(hObject, eventdata, handles)
% hObject    handle to sig_thresh_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fc_static_menu_Callback(hObject, eventdata, handles)
% hObject    handle to fc_static_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fc_dyn_menu_Callback(hObject, eventdata, handles)
% hObject    handle to fc_dyn_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function surrogate_menu_Callback(hObject, eventdata, handles)
% hObject    handle to surrogate_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function permutation_menu_Callback(hObject, eventdata, handles)
% hObject    handle to permutation_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function pca_dim_reduc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to pca_dim_reduc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fc_qc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to fc_qc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function edge_group_glm_menu_Callback(hObject, eventdata, handles)
% hObject    handle to edge_group_glm_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function data_qc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to data_qc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function util_menu_Callback(hObject, eventdata, handles)
% hObject    handle to util_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function gen_template_menu_Callback(hObject, eventdata, handles)
% hObject    handle to gen_template_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function settings_menu_Callback(hObject, eventdata, handles)
% hObject    handle to settings_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function dcp_import_menu_Callback(hObject, eventdata, handles)
% hObject    handle to dcp_import_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function add_aux_data_menu_Callback(hObject, eventdata, handles)
% hObject    handle to add_aux_data_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_subj_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_subj_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_results_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_results_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function net_group_glm_menu_Callback(hObject, eventdata, handles)
% hObject    handle to net_group_glm_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function create_group_proj_menu_Callback(hObject, eventdata, handles)
% hObject    handle to create_group_proj_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function physio_menu_Callback(hObject, eventdata, handles)
% hObject    handle to physio_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global X Y Z sagslice Vflip
sagslice=get(handles.slider1,'Value');cla(handles.axes3)
if ~isempty(Vflip)
    contourslice(handles.axes3,Y,X,Z,Vflip,[],sagslice,[],3);
    set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
end

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global X Y Z axslice Vflip
axslice=get(handles.slider2,'Value');cla(handles.axes4)
if ~isempty(Vflip)
    contourslice(handles.axes4,Y,X,Z,Vflip,[],[],axslice,3);
    set(handles.axes4,'view',[-90 90],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
end


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global X Y Z corslice Vflip
corslice=get(handles.slider3,'Value');cla(handles.axes5)
if ~isempty(Vflip)
    contourslice(handles.axes5,Y,X,Z,Vflip,corslice,[],[],3);
    set(handles.axes5,'view',[-90 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
end


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_fc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_fc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_data dcp_data_orig circ_thresh isdfc dfc_td

[thisfile,~]=uigetfile; %('*.mat');
if thisfile~=0
    [~,~,fext]=fileparts(thisfile);
end
if thisfile~=0 && strcmp(fext,'.txt')
    dcp_data_orig=load(thisfile);
    dcp_data=dcp_data_orig;
end
if thisfile~=0 && strcmp(fext,'.mat')
    varname=cell2mat(who('-file',thisfile));
    load(thisfile,varname);
    dcp_data_orig=varname;
    dcp_data=dcp_data_orig;
end

if thisfile~=0
    fcsize=size(dcp_data);
    if ndims(fcsize)>1
        isdfc=1;
        dfc_td=fcsize(1);
    else
        isdfc=0;
    end

    if isdfc~=1
        imagesc(handles.axes1,squeeze(dcp_data(:,:)));colorbar;colormap('jet')
        % set(gca,'YTick',0:length_data-1);
        % yticklabels(dcp_nodelabels);
        % set(gca,'FontSize',3);ytickangle(35);
        dcp_data_graph=squeeze(dcp_data(1,:,:));
        dcp_data_graph(dcp_data_graph<circ_thresh)=0;
        G=graph(dcp_data_graph,'upper');
        axes(handles.axes2);cla
        plot(handles.axes2,G,'Layout','circle')
    end

    if isdfc==1
        imagesc(handles.axes1,squeeze(dcp_data(1,:,:)));colorbar;colormap('jet')
        % set(gca,'YTick',0:length_data-1);
        % yticklabels(dcp_nodelabels);
        % set(gca,'FontSize',3);ytickangle(35);
        dcp_data_graph=squeeze(dcp_data(1,:,:));
        dcp_data_graph(dcp_data_graph<circ_thresh)=0;
        G=graph(dcp_data_graph,'upper');
        axes(handles.axes2);cla
        plot(handles.axes2,G,'Layout','circle')
    end
end


% --------------------------------------------------------------------
function load_dfc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_dfc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_data dcp_data_orig
[thisfile,~]=uigetfile; %('*.mat');
if thisfile~=0
    dcp_data_orig=load(thisfile);
end


% --------------------------------------------------------------------
function load_coords_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_coords_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_coords dcp_coords_vox dcp_coords_x dcp_coords_y dcp_coords_z
[thisfile,~]=uigetfile; %('*.mat');
if thisfile~=0
    coordvar=cell2mat(who('-file',thisfile));
    load(thisfile);
    dcp_coords=eval(coordvar);
    dcp_coords_vox(:,2)=dcp_coords(:,2)+127; % Y
    dcp_coords_vox(:,3)=dcp_coords(:,3)+79; % Z
    dcp_coords_vox(:,1)=dcp_coords(:,1)+90; % X
    dcp_coords_x=mat2cell(dcp_coords(:,1),ones(1,size(dcp_coords,1)));
    dcp_coords_y=mat2cell(dcp_coords(:,2),ones(1,size(dcp_coords,1)));
    dcp_coords_z=mat2cell(dcp_coords(:,3),ones(1,size(dcp_coords,1)));
end


% --------------------------------------------------------------------
function load_labels_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_labels_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_nodelabels_cell dcp_coords_x dcp_coords_y dcp_coords_z dcp_coords
[thisfile,~]=uigetfile; %('*.mat');
if thisfile~=0
    labelvar=cell2mat(who('-file',thisfile));
    load(thisfile);
    dcp_nodelabels_cell=eval(labelvar);
    if ~isempty(dcp_coords)
        set(handles.uitable1,'Data',[dcp_nodelabels_cell dcp_coords_x dcp_coords_y dcp_coords_z]);
    else
        set(handles.uitable1,'Data',[dcp_nodelabels_cell]);
    end
end


% --------------------------------------------------------------------
function load_braintemplate_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_braintemplate_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_template V X Y Z sagslice axslice corslice Vflip
% [thisfile,~]=uigetfile; %('*.mat');
% if thisfile~=0
%     dcp_template=load(thisfile);
% end

    [~,V,~,~]=BrikLoad('MNI152_T1_2mm+tlrc');
    [Y,X,Z]=meshgrid(1:109,1:91,1:91);
    
    Vflip=V; %flip(V,2);
    
    sagslice=45;axslice=45;corslice=54;
    
%     dcp_coords_vox(:,2)=dcp_coords(:,2)+127; % Y
%     dcp_coords_vox(:,3)=dcp_coords(:,3)+79; % Z
%     dcp_coords_vox(:,1)=dcp_coords(:,1)+90; % X

    % Generate Sagittal slice images
    contourslice(handles.axes3,Y,X,Z,Vflip,[],45,[],3)
    set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
    daspect(handles.axes3,[1 1 1]);
    
    % Generate Axial slice images
    contourslice(handles.axes4,Y,X,Z,Vflip,[],[],45,3)
    set(handles.axes4,'view',[-90 90],'XLim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
    daspect(handles.axes4,[1 1 1]);
    
    % Generate Coronal slice images
    contourslice(handles.axes5,Y,X,Z,Vflip,54,[],[],3)
    set(handles.axes5,'view',[-90 0],'XLim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
    daspect(handles.axes5,[1 1 1]);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --------------------------------------------------------------------
function load_net_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_net_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_nodelabels_cell dcp_icnid dcp_icnlabels
[thisfile,~]=uigetfile; %('*.mat');
if thisfile~=0
    dcp_icnlabels=load(thisfile);
    dcp_icnid=load(thisfile);
end


% --------------------------------------------------------------------
function get_fc_menu_Callback(hObject, eventdata, handles)
% hObject    handle to get_fc_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% keyboard
% wsvars=evalin('base','who');
% if isempty(wsvars)
%     wsvars='none';
% end
% handles.popupmenu2.String=wsvars;


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
global dcp_data dcp_data_orig


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_data_orig dcp_data isdfc dfc_td circ_thresh tframe Vflip dcp_coords_vox
global X Y Z sagslice axslice corslice numrois

circ_thresh=str2num(handles.thresh_editbox.String);
if isempty(tframe)
    tframe=1;
end

if ~strcmp(handles.popupmenu2.String,'none')
    varname=handles.popupmenu2.String{handles.popupmenu2.Value};
    dcp_data_orig=evalin('base',varname);
    dcp_data=dcp_data_orig;

    fcsize=size(dcp_data);
    numrois=fcsize(2);
    if ndims(dcp_data)>2
        isdfc=1;
        dfc_td=fcsize(1);
    else
        isdfc=0;
    end

    if isdfc~=1
        circ_thresh=.75;set(handles.thresh_editbox,'String',num2str(circ_thresh));
        imagesc(handles.axes1,dcp_data);colorbar(handles.axes1);colormap('jet')
        % set(gca,'YTick',0:length_data-1);
        % yticklabels(dcp_nodelabels);
        % set(gca,'FontSize',3);ytickangle(35);
        dcp_data_graph=dcp_data;
        dcp_data_graph(dcp_data_graph<circ_thresh)=0;
        G=graph(dcp_data_graph,'upper');
        axes(handles.axes2);cla
        plot(handles.axes2,G,'Layout','circle')
        tframe=1;
        set(handles.time_editbox,'String',tframe);
        
%         dcp_data_graph(dcp_data_graph==tril(dcp_data_graph))=0;
        [I,J]=find(dcp_data_graph>=circ_thresh);
        if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
            FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
            cla(handles.axes3)
            contourslice(handles.axes3,Y,X,Z,Vflip,[],sagslice,[],3);
            set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
            hold(handles.axes3,'on')
            for loop1=1:size(FF,1)
                plot3(handles.axes3,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
            end
    %     daspect(handles.axes3,[1 1 1]);
        end
        
        if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
            [I,J]=find(dcp_data_graph>=circ_thresh);
            FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
            cla(handles.axes4)
            contourslice(handles.axes4,Y,X,Z,Vflip,[],[],axslice,3);
            set(handles.axes4,'view',[-90 90],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
            hold(handles.axes4,'on')
            for loop1=1:size(FF,1)
                plot3(handles.axes4,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
            end
        %     daspect(handles.axes4,[1 1 1]);
        end

        if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
            [I,J]=find(dcp_data_graph>=circ_thresh);
            FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
            cla(handles.axes5)
            contourslice(handles.axes5,Y,X,Z,Vflip,corslice,[],[],3);
            set(handles.axes5,'view',[-90 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
            hold(handles.axes5,'on')
            for loop1=1:size(FF,1)
                plot3(handles.axes5,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
            end
        %     daspect(handles.axes5,[1 1 1]);
        end
    end
    
    if isdfc==1
        circ_thresh=.75;set(handles.thresh_editbox,'String',num2str(circ_thresh));
        imagesc(handles.axes1,squeeze(dcp_data(tframe,:,:)));colorbar(handles.axes1);colormap('jet')
        % set(gca,'YTick',0:length_data-1);
        % yticklabels(dcp_nodelabels);
        % set(gca,'FontSize',3);ytickangle(35);
        dcp_data_graph=squeeze(dcp_data(tframe,:,:));
        dcp_data_graph(dcp_data_graph<circ_thresh)=0;
        G=graph(dcp_data_graph,'upper');
        axes(handles.axes2);cla
        plot(handles.axes2,G,'Layout','circle')
        set(handles.time_editbox,'String',handles.time_editbox.String);
        
        dcp_data_graph(dcp_data_graph==tril(dcp_data_graph))=0;
        [I,J]=find(dcp_data_graph>=circ_thresh);
        if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
            FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
            cla(handles.axes3)
            contourslice(handles.axes3,Y,X,Z,Vflip,[],sagslice,[],3);
            set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
            hold(handles.axes3,'on')
            for loop1=1:size(FF,1)
                plot3(handles.axes3,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
            end
    %     daspect(handles.axes3,[1 1 1]);
        end
        
        if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
            [I,J]=find(dcp_data_graph>=circ_thresh);
            FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
            cla(handles.axes4)
            contourslice(handles.axes4,Y,X,Z,Vflip,[],[],axslice,3);
            set(handles.axes4,'view',[-90 90],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
            hold(handles.axes4,'on')
            for loop1=1:size(FF,1)
                plot3(handles.axes4,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
            end
        %     daspect(handles.axes4,[1 1 1]);
        end

        if ~isempty(dcp_coords_vox) && ~isempty(Vflip)
            [I,J]=find(dcp_data_graph>=circ_thresh);
            FF=dcp_coords_vox(I,:)./2;TT=dcp_coords_vox(J,:)./2;
            cla(handles.axes5)
            contourslice(handles.axes5,Y,X,Z,Vflip,corslice,[],[],3);
            set(handles.axes5,'view',[-90 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
            hold(handles.axes5,'on')
            for loop1=1:size(FF,1)
                plot3(handles.axes5,[FF(loop1,2);TT(loop1,2)],[FF(loop1,1);TT(loop1,1)],[FF(loop1,3);TT(loop1,3)],'b','LineWidth',1.5,'Marker','o','MarkerFaceColor','b')
            end
        %     daspect(handles.axes5,[1 1 1]);
        end
    end
    
end


% --------------------------------------------------------------------
function utils_menu_Callback(hObject, eventdata, handles)
% hObject    handle to utils_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function make_template_menu_Callback(hObject, eventdata, handles)
% hObject    handle to make_template_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of togglebutton1
global axslice corslice sagslice dcp_data circ_thresh Vflip X Y Z
% keyboard
if handles.togglebutton1.Value==0
    handles.togglebutton1.String='Anat.';
    imagesc(handles.axes3,squeeze(Vflip(sagslice,:,:)));colormap('gray');view(handles.axes3,[-90 90])
    imagesc(handles.axes4,squeeze(Vflip(:,:,axslice)));colormap('gray');view(handles.axes4,[-90 90])
    imagesc(handles.axes5,squeeze(Vflip(:,corslice,:)));colormap('gray');view(handles.axes5,[90 -90])
end
if handles.togglebutton1.Value==1
    handles.togglebutton1.String='Contour';
    
    % Generate Sagittal slice images
    cla(handles.axes3)
    contourslice(handles.axes3,Y,X,Z,Vflip,[],sagslice,[],3);colormap('parula')
    set(handles.axes3,'view',[0 0],'Xlim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
    daspect(handles.axes3,[1 1 1]);
    
    % Generate Axial slice images
    cla(handles.axes4)
    contourslice(handles.axes4,Y,X,Z,Vflip,[],[],axslice,3);colormap('parula')
    set(handles.axes4,'view',[-90 90],'XLim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
    daspect(handles.axes4,[1 1 1]);
    
    % Generate Coronal slice images
    cla(handles.axes5)
    contourslice(handles.axes5,Y,X,Z,Vflip,corslice,[],[],3);colormap('parula')
    set(handles.axes5,'view',[-90 0],'XLim',[1 109],'YLim',[1 91],'ZLim',[1 81],'YTick',[],'XTick',[],'ZTick',[]);
    daspect(handles.axes5,[1 1 1]);
end
% keyboard

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
wsvars=evalin('base','who');
if isempty(wsvars)
    wsvars='none';
end
handles.popupmenu2.String=wsvars;



function keeproieditbox_Callback(hObject, eventdata, handles)
% hObject    handle to keeproieditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of keeproieditbox as text
%        str2double(get(hObject,'String')) returns contents of keeproieditbox as a double


% --- Executes during object creation, after setting all properties.
function keeproieditbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to keeproieditbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_data dcp_data_orig dcp_nodelabels_cell numrois isdfc

if ~strcmp(handles.keeproieditbox.String,'Keep ROIs')
    if ~isempty(dcp_data) % && ~isempty(dcp_nodelabels_cell)
        roi_list=str2num(handles.keeproieditbox.String);
        if ~isempty(roi_list)
            rmroi=1:1:numrois;
            for loop1=1:length(roi_list)
                rmroi(rmroi==roi_list(loop1))=[];
            end
            if isdfc==1
                dcp_data(:,rmroi,rmroi)=0;
            else
                dcp_data(rmroi,rmroi)=0;
            end
        end
    end
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dcp_nodelabels_cell dcp_data circ_thresh nodewt tframe dcp_coords_x dcp_coords_y dcp_coords_z

if ~isempty(dcp_nodelabels_cell) && ~isempty(dcp_data)
    dcp_data_thr=dcp_data;
    dcp_data_thr(dcp_data_thr<circ_thresh)=0;
    nodewt=cellstr(num2str(sum(squeeze(dcp_data_thr(tframe,:,:)),2)));
    set(handles.uitable1,'Data',[dcp_nodelabels_cell dcp_coords_x dcp_coords_y dcp_coords_z cell(size(dcp_coords_x)) nodewt]);
end


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% % % % % % % % % % % % % 
