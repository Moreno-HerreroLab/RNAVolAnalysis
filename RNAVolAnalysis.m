function varargout = RNAVolAnalysis(varargin)
% RNAVOLANALYSIS MATLAB code for RNAVolAnalysis.fig
%      RNAVOLANALYSIS, by itself, creates a new RNAVOLANALYSIS or raises the existing
%      singleton*.
%
%      H = RNAVOLANALYSIS returns the handle to a new RNAVOLANALYSIS or the handle to
%      the existing singleton*.
%
%      RNAVOLANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RNAVOLANALYSIS.M with the given input arguments.
%
%      RNAVOLANALYSIS('Property','Value',...) creates a new RNAVOLANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RNAVolAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RNAVolAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help RNAVolAnalysis
% Last Modified by GUIDE v2.5 05-Jun-2023 15:14:10
% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RNAVolAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @RNAVolAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT



% --- Executes during object creation, after setting all properties.
function data_path_CreateFcn(hObject, ~, ~)

%data_directory = 'C:\Users\evatr\Documents\Script_for_rna_analysis\datasets\SARS_35';
data_directory = 'C:\Users\Eva Martin\Documents\MorenoHerreroLab\Projects_tesis\SARS\Analisis SARS para paper\IMAGING DATA\Sars-CoV-2\individual molecules';
set(hObject,'String', num2str(data_directory));

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes just before RNAVolAnalysis is made visible.
function RNAVolAnalysis_OpeningFcn(hObject, ~, handles, varargin)

%Here we inizialite the values of the GUI
set(handles.smooth_edit, 'String', '5')
set(handles.sort_edit,'Max', 50)
set(handles.sort_edit,'String', 'Sorted List')
handles.Results_Cell={};
handles.color='b';

handles.resize_factor=512;
handles.size_nm=80;

set(handles.pixel_size,'String', num2str(handles.resize_factor));
set(handles.nm_size,'String', num2str(handles.size_nm));

handles.factor=handles.size_nm/handles.resize_factor;

% Choose default command line output for RNAVolAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RNAVolAnalysis wait for user response (see UIRESUME)
%uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RNAVolAnalysis_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browser.
function browser_Callback(~, ~, handles)

input_pathname=uigetdir();
set(handles.data_path,'String', num2str(input_pathname));



function data_path_Callback(hObject, ~, ~)

chckfldr=get(hObject,'String');
chckfldr= exist(chckfldr,'dir');
if chckfldr ~= 7,  msgbox('The provided directory is not valid. Please browse or insert a valid one','ERROR', 'error')
     return; end


% --- Executes on button press in load_first.
function load_first_Callback(hObject, ~, handles)

handles.order=0;
handles.data_directory=get(handles.data_path,'String');
CurrentFolder=pwd;
handles.fileindex=1;
cd(handles.data_directory);

handles.AllFileNames=dir('*.txt');
handles.JPGFileNames=dir('*.jpg');
cd(CurrentFolder)
handles.chckflder=length(handles.AllFileNames);

if handles.chckflder == 0
    msgbox('The provided input folder is empty.','ERROR', 'error')
else 
    [handles.RNA, handles.FileName]= func_load(handles.data_directory,handles.AllFileNames, handles.fileindex,handles.resize_factor);
    isempty(handles.JPGFileNames);
    handles.RNAjpg=[];
    if isempty(handles.JPGFileNames) == 0
        [handles.RNAjpg, handles.JPGFileName]= func_load_jpg(handles.data_directory,handles.JPGFileNames, handles.fileindex,handles.resize_factor);
    end 
    set(handles.frame_name,'String', num2str(handles.FileName));       
    axes(handles.axes1);
    cla reset
    imshow(handles.RNA);
        
    max_val=max(max(handles.RNA))+0.0001;
    min_val=min(min(handles.RNA))-0.0001;
    
    set(handles.thres_slider,'Min',min_val,'Max',max_val,'SliderStep',[0.01, 0.1],'Value',min_val+0.1)
    set(handles.threshold_edit,'String',num2str(min_val+0.1));

end 

guidata(hObject, handles);


% --- Executes on button press in load_next.
function load_next_Callback(hObject, ~, handles)

handles.order=0;
if handles.fileindex == handles.chckflder
    msgbox('No next frame','ERROR', 'error')
else
    handles.fileindex=handles.fileindex+1;
    [handles.RNA, handles.FileName]= func_load(handles.data_directory,handles.AllFileNames, handles.fileindex,handles.resize_factor);
    handles.RNAjpg=[];
    set(handles.frame_name,'String', num2str(handles.FileName));
    if isempty(handles.JPGFileNames) == false
        [handles.RNAjpg, handles.JPGFileName]= func_load_jpg(handles.data_directory,handles.JPGFileNames, handles.fileindex,handles.resize_factor);

    end 

    axes(handles.axes1);
    cla reset
    imshow(handles.RNA);
    
    max_val=max(max(handles.RNA))+0.0001;
    min_val=min(min(handles.RNA))-0.0001;

    
    set(handles.thres_slider,'Min',min_val,'Max',max_val,'SliderStep',[0.01, 0.1],'Value',min_val+0.1);
    set(handles.threshold_edit,'String',num2str(min_val+0.1));


end

guidata(hObject, handles);

% --- Executes on button press in load_prev.
function load_prev_Callback(hObject, ~, handles)

handles.order=0;
if handles.fileindex == 1
    msgbox('No pevious frame','ERROR', 'error')
else
    handles.fileindex=handles.fileindex-1;
    [handles.RNA, handles.FileName]= func_load(handles.data_directory,handles.AllFileNames, handles.fileindex,handles.resize_factor);
    handles.RNAjpg=[];
    set(handles.frame_name,'String', num2str(handles.FileName));
    if isempty(handles.JPGFileNames) == false
        [handles.RNAjpg, handles.JPGFileName]= func_load_jpg(handles.data_directory,handles.JPGFileNames, handles.fileindex,handles.resize_factor);

    end 

    axes(handles.axes1);
    cla reset
    imshow(handles.RNA);

    max_val=max(max(handles.RNA))+0.0001;
    min_val=min(min(handles.RNA))-0.0001;
  
    set(handles.thres_slider,'Min',min_val,'Max',max_val,'SliderStep',[0.01, 0.1],'Value',min_val+0.1);
    set(handles.threshold_edit,'String',num2str(min_val+0.1));



end

guidata(hObject, handles);

% --- Executes on button press in smooth_button.
function smooth_button_Callback(hObject, ~, handles)

smooth_param=str2double(get(handles.smooth_edit,'String'));
handles.RNA2=imgaussfilt(handles.RNA,smooth_param);

axes(handles.axes2);
cla reset 
imshow(handles.RNA2)

guidata(hObject, handles);



function smooth_edit_CreateFcn(hObject, ~, ~)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function smooth_edit_Callback(hObject, ~, handles)


guidata(hObject, handles);

% --- Executes on slider movement.
function thres_slider_Callback(hObject, ~, handles)

handles.threshold=get(handles.thres_slider,'Value');
set(handles.threshold_edit,'String',num2str(handles.threshold));
RNA3=apply_threshold(handles);
handles.RNA3=RNA3;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function thres_slider_CreateFcn(hObject, ~, ~)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in push_auto_thres.
function push_auto_thres_Callback(hObject, ~, handles)

handles.threshold=isodata(handles.RNA2);
set(handles.threshold_edit,'String',num2str(handles.threshold));
set(handles.thres_slider,'Value',handles.threshold);

handles.RNA3=apply_threshold(handles);

guidata(hObject, handles);


% --- Executes on button press in run.
function run_Callback(hObject, ~, handles)

[linkers,objects]=main_run(handles);
handles.linkers=linkers;
handles.objects=objects;
[handles.Volume,handles.Noise,handles.labelMatrix]=extend_molecule(handles);

cell1={handles.Volume.vol};
cell2={handles.Volume.info};
handles.TAB=[cell1;cell2];

set(handles.uitable1,'Data',handles.TAB')

handles.RNA_profile=zeros(size(handles.RNA2));


if isempty(handles.RNAjpg) == 0 
    axes(handles.axes4)
    cla reset
    imshow(handles.RNAjpg)

else
    axes(handles.axes4)
    cla reset
    imshow(handles.RNA)

end 

guidata(hObject, handles);


% --- Executes on button press in select_thres.
function select_thres_Callback(hObject, ~, handles)

handles.mol_thres=handles.threshold;
handles.thres_array=[handles.mol_thres];
handles.vol_sort_array=[];
handles.info_sort_array={};


guidata(hObject, handles);

% --- Executes on button press in select_domains_thres.
function select_domains_thres_Callback(hObject, ~, handles)

thres_n=handles.threshold;
handles.thres_array=[handles.thres_array,thres_n];

guidata(hObject, handles);


% --- Executes on button press in sort_button.
function sort_button_Callback(hObject, ~, handles)

handles.order=handles.order+1;
vol=handles.Volume(handles.row).vol;
info=handles.Volume(handles.row).info;



handles.vol_sort_array=[handles.vol_sort_array, vol];
handles.info_sort_array{handles.order}=info;

s=string(handles.vol_sort_array);
set(handles.sort_edit,'String',s)

handles.data=str2double(s);
handles.thres=[handles.thres_array(1) handles.thres_array(end)];
handles.infor = string(handles.info_sort_array);

axes(handles.axes4);
hold on
[x,y]=find(handles.labelMatrix == handles.row);
%s=scatter(y,x,'filled','SizeData',1);
s=scatter(y,x,'filled','SizeData',1,'MarkerFaceColor',handles.color);
s.MarkerFaceAlpha = .3;

guidata(hObject, handles);


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)

IndicesCopia=eventdata.Indices;
handles.row=IndicesCopia(:,1); %double 
%handles.row_info=IndicesCopia(:,2); %double 
guidata(hObject, handles);


% --- Executes on button press in clear_sort_button.
function clear_sort_button_Callback(hObject, ~, handles)

handles.vol_sort_array=[];
handles.info_sort_array={};
handles.order=0;
s=[];
set(handles.sort_edit,'String',s);
handles.RNA_profile=zeros(size(handles.RNA2));

guidata(hObject,handles)

% --- Executes on selection change in popupmenu_colores.
function popupmenu_colores_Callback(hObject, ~, handles)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_colores contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_colores

C=get(hObject,'Value');
if C==1
    handles.color='b';
elseif C==2
    handles.color='g';
elseif C==3
    handles.color='r';
elseif C==4
    handles.color='c';
elseif C==5
    handles.color='m';
elseif C==6
    handles.color='y';
elseif C==7
    handles.color='k';
elseif C==8
    handles.color='w';
end 
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu_colores_CreateFcn(hObject, ~, ~)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function region_color_edit_Callback(~, ~, ~)

% Hints: get(hObject,'String') returns contents of region_color_edit as text
%        str2double(get(hObject,'String')) returns contents of region_color_edit as a double


% --- Executes during object creation, after setting all properties.
function region_color_edit_CreateFcn(hObject, ~, ~)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in change_color_button.
function change_color_button_Callback(hObject, ~, handles)
% hObject    handle to change_color_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.reg=str2double(get(handles.region_color_edit,'String'));
axes(handles.axes4);
hold on

[x,y]=find(handles.labelMatrix == handles.reg);

s=scatter(y,x,'filled','SizeData',1,'MarkerFaceColor',handles.color);
s.MarkerFaceAlpha = .2;


% label=[label;num2str(handles.reg)];
% legend(label)
% legend

guidata(hObject, handles);




% --- Executes on button press in save_molecule.
function save_molecule_Callback(hObject, ~, handles)
% hObject    handle to save_molecule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


st=handles.FileName;
ix=find(st=='.');
num=str2double(st(1:ix-1));


cell1={num, handles.vol_sort_array};
cellinfo={'info', handles.infor};
cell2={'thresholds',handles.thres_array};
cell3={'background noise',handles.Noise};

handles.Results_Cell=[handles.Results_Cell;cell1,cellinfo,cell2,cell3];


writecell(handles.Results_Cell,'Results.dat');

guidata(hObject, handles);



function sort_edit_Callback(~, ~, ~)

% Hints: get(hObject,'String') returns contents of sort_edit as text
%        str2double(get(hObject,'String')) returns contents of sort_edit as a double


% --- Executes during object creation, after setting all properties.
function sort_edit_CreateFcn(hObject, ~, ~)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold_edit_Callback(hObject, ~, handles)

handles.threshold=str2double(get(handles.threshold_edit,'String'));
set(handles.thres_slider,'Value',handles.threshold)
RNA3=apply_threshold(handles);
handles.RNA3=RNA3;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function threshold_edit_CreateFcn(hObject, ~, ~)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function nm_size_Callback(hObject, eventdata, handles)
% hObject    handle to nm_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nm_size as text
%        str2double(get(hObject,'String')) returns contents of nm_size as a double


handles.size_nm=str2double(get(handles.nm_size,'String'));
handles.factor=handles.size_nm/handles.resize_factor;


guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nm_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nm_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pixel_size_Callback(hObject, eventdata, handles)
% hObject    handle to pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_size as text
%        str2double(get(hObject,'String')) returns contents of pixel_size as a double

handles.resize_factor=str2double(get(handles.pixel_size,'String'));
handles.factor=handles.size_nm/handles.resize_factor;


guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function pixel_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function RNA1=read_RemoveHeader(path)
fid1 = fopen(path,'r');
linedata = '0';
lineNumber = 0;
mat=[];

while (ischar(linedata))
    linedata = fgetl(fid1);
    lineNumber = lineNumber +1;
    if lineNumber > 5
        if linedata~=-1
        x = str2num(linedata);
        mat=[mat;x];
        end 
        
    end
end

fclose(fid1);
newImage = flip(mat,1);
RNA1 = flip(newImage,2);

function RNA3=apply_threshold(handles)

RNA3=handles.RNA2;
mask=handles.RNA2>handles.threshold;
RNA3(mask==0)=0;
axes(handles.axes2)
cla reset
imshow(RNA3);

function [linkers,objects]=main_run(handles)

handles.thres_array=sort(handles.thres_array);

n_steps=length(handles.thres_array);
linker_map=zeros(size(handles.RNA2));
object_map=zeros(size(handles.RNA2));

for i=1:n_steps-1
    low=handles.thres_array(i);
    high=handles.thres_array(i+1);
    
    %guardamos la lista de pixeles de cada isla con el thres alto
    maskhigh=zeros(size(handles.RNA2));
    maskhigh(handles.RNA2 > high) = 1; 
    maskhigh=imclearborder(maskhigh);
    CC_high = bwconncomp(maskhigh);    
    idx_list_high = CC_high.PixelIdxList;
    
    %guardamos la lista de pixeles de cada isla con el thres bajo     
    masklow=zeros(size(handles.RNA2));
    masklow(handles.RNA2 > low) = 1;
    masklow=imclearborder(masklow);
    CC_low = bwconncomp(masklow);   
    idx_list_low = CC_low.PixelIdxList;  

    %comparamos cada lista de pixeles bajo con cada lista de pixeles en el thres
    %alto. Es decir, cada isla de thres bajo tiene de opciones:
    
    %A. no coincidir con ninguna isla del alto
    %B. coincidir con una isla del alto
    %C. coincidir con 2 o más islas del alto
    
    %En los casos A y B no se hace nada,  en el caso C: guardamos la lista
    %de pixeles de las regiones que van a unirse (thres alto) y vamos
    %bajando el threshold hasta que se unan. Una vez unidas, se hace el
    %skeleton. Al skeleton se le restan los dominios de thres_alto, y nos
    %quedamos con el linker que toca las regiones.
    
    for k=1:length(idx_list_low)        
        object_n=idx_list_low{k};
        tf = [];
        for j=1:length(idx_list_high)
            object_n2=idx_list_high{j};
            C = ismember(object_n,object_n2); %C va a ser un array con 1 en las posiciones que coincidan y 0 en las que no
            %Si C es todo ceros: el objeto n no coincide con el n2
            %Si C esta compuesto de ceros y unos: el objeto si coincide
            is_new=find(C==1, 1); %Este vector estará vacío si las islas no coinciden, y no lo estará si las isla coinciden 
            if isempty(is_new)==1
                tf = [tf, 0]; %0 indica que el objeto no coincide 
            else
                tf=[tf, 1]; %1 indica que el objeto coincide 
            end 
        end 
        %Si el tf hay mas de un 1: Caso C
        %Si en el tf no hay ningún 1: caso A
        %Si en el tf solo hay un 1: Caso B
        val=sum(tf);
        if val >1
%            %guardamos las islas que se han unido
             prev_objects_idx=find(tf==1);
             
%            %objetos se han unido en uno solo. SACAR LINKER 
             for  indexes=1:length(prev_objects_idx)
             object_unido_n_idxs=idx_list_high{prev_objects_idx(indexes)};
             object_map_2=zeros(size(handles.RNA2));
             object_map_2(object_unido_n_idxs)=1;
             
             object_map=object_map + object_map_2;
             
             end 

            link= zeros(size(handles.RNA2));
            link(object_n)=1;
            skel = bwmorph(link,'thin',Inf);
            
            whole=maskhigh+skel;
            whole(whole>0)=1;
            linkers=whole-maskhigh;

            %Nos quedamos solo con los branches importantes
            %----------------------------------------------------------------------------------
            CLink=bwconncomp(linkers);        
            linkers_def=zeros(size(handles.RNA2));

            for h=1:CLink.NumObjects
                count=0;
                link_n=zeros(size(handles.RNA2));
                link_n(CLink.PixelIdxList{h})=1;
                ep=bwmorph(link_n,'endpoints');
                ep_loc=find(ep);
                for ep_i=1:length(ep_loc)
                    [~,mask]=find_8_conn(ep_loc(ep_i),handles.RNA2);
                    aux=mask-maskhigh;
                    sur=find(aux==1);
                    if length(sur)<8 %este ep toca alguna region
                        count=count+1;
                    end
                end
                if count>1
                    linkers_def=linkers_def+link_n;
                end
            end 

         link_domain=extend_border_link(linkers_def, high-0.01 , low , handles.RNA2);
         linker_map=linker_map+link_domain;
        end
    end

end 

L= max(object_map,[],'all');
object_map_def=zeros(size(handles.RNA2));

for i=1:L
    aux=zeros(size(handles.RNA2));
    [ind_rep]=find(object_map == i);
    aux(ind_rep)=1;
    aux = bwpropfilt(logical(aux),'EulerNumber',[1 1]);
    object_map_def=object_map_def+aux;
    
end 

axes(handles.axes2);
cla reset
imshow(handles.RNA)
hold on 
whole_border=bwperim(object_map_def);
[x, y] = find(whole_border);

plot(y,x,'.','LineWidth',20);

[x,y]=find(linker_map);
plot(y,x,'r.','LineWidth',20);

linkers=linker_map;
objects=object_map_def;




function [Volume,basal_plane_height,labelMatrix]=extend_molecule(handles)

linkers=handles.linkers;
CCL=bwconncomp(linkers);
objects=handles.objects;
CCO=bwconncomp(objects);
lista_linkers=CCL.PixelIdxList;
lista_objects=CCO.PixelIdxList;
total_dom=length(lista_objects) + length(lista_linkers);

handles.labelMatrix=zeros(size(handles.RNA2));


for i=1:length(lista_objects)    
    handles.labelMatrix(lista_objects{i})=i;
end 
for i=length(lista_objects)+1:total_dom    
    indx= i-(length(lista_objects));
    handles.labelMatrix(lista_linkers{indx})=i;

end

for g=1:200
whole_border=bwperim(handles.labelMatrix);
positions=find(whole_border);
[row, colum] = find(whole_border);
[r_total,c_total]=size(handles.RNA2);
RNA2=handles.RNA2;
for i=1:length(positions)
    pos=positions(i);
    int = handles.labelMatrix(pos);
    c=colum(i);
    r=row(i);
    if r+1 < r_total && r-1 > 0 && c+1 < c_total && c-1 >0
    izq = handles.labelMatrix(r, c +1);
    der = handles.labelMatrix(r, c-1);
    arr = handles.labelMatrix(r+1,c);
    abj = handles.labelMatrix(r-1,c);
    
    
    if izq ==0 && RNA2(r, c +1) > handles.mol_thres %The mask condition is to stop when the molecule ends
        handles.labelMatrix(r, c +1) = int;
    end
    if der ==0 && RNA2(r, c +1) > handles.mol_thres
        handles.labelMatrix(r, c -1) = int;
    end
    if arr ==0 && RNA2(r, c +1) >handles.mol_thres
        handles.labelMatrix(r+1, c) = int;
    end
    if abj ==0 && RNA2(r, c +1) > handles.mol_thres
        handles.labelMatrix(r-1, c) = int;
    end
    end 

end 
end 

axes(handles.axes3);
hold on 
handles.RNA_fig=handles.RNA;
%label=[];
for i=1:total_dom
    [x,y]=find(handles.labelMatrix == i);   
    handles.RNA_fig=insertText(handles.RNA_fig, [mean(y) mean(x)], num2str(i),'FontSize',20);
end 
imshow(handles.RNA_fig)

for i=1:total_dom
    positt=handles.labelMatrix == i;
    aux=zeros(size(handles.RNA));
    aux(positt)=1;
    border=bwperim(aux);
    [xx,yy]=find(border);
    plot(yy,xx,'.r','MarkerSize',3);
    
end 



noise_locs=find(handles.RNA2<handles.mol_thres);
noise_values=handles.RNA2(noise_locs);
basal_plane_height=mean(noise_values);
handles.basal_plane=basal_plane_height;
total_volume=0;
Volume=[];
volume_noise=0;
for i=1:total_dom
    ind=find(handles.labelMatrix == i);
    values=handles.RNA2(ind)-basal_plane_height;
    
    values_noise=handles.RNA2(ind)-basal_plane_height;

    volume=sum(values)*((handles.factor)^2); %factor 
    volume_noise=volume_noise + sum(values_noise)*((handles.factor)^2); %factor

    total_volume=total_volume + volume;
    Volume(i).vol=volume;

    if i<= length(lista_objects)
        Volume(i).info='Domain';
    else
        Volume(i).info='Linker';
    end 
end 


labelMatrix=handles.labelMatrix;




function [extended_bp_locations,mask]=find_8_conn(b_locs,RNA)

mask=zeros(size(RNA));
mask(b_locs)=1;

[r,c]=find(mask);
mask=zeros(size(RNA));

    mask(r+1,c)=1;
    mask(r-1,c)=1;
    mask(r,c+1)=1;
    mask(r,c-1)=1;
    mask(r-1,c-1)=1;
    mask(r+1,c+1)=1;
    mask(r+1,c-1)=1;
    mask(r-1,c+1)=1;
    
mask=imclearborder(mask);
extended_bp_locations=find(mask);

function [link_domain2]=extend_border_link(link,thres_alto,thres_bajo,RNA2)

[r_total,c_total]=size(RNA2);

for h=1:1
whole_border=bwperim(link);
[row, colum] = find(whole_border);
positions = find(whole_border);
    for i=1:length(positions)
        c=colum(i);
        r=row(i);
        if r+1 < r_total && r-1 > 0 && c+1 < c_total && c-1 >0
            if RNA2(r+1,c+1) < thres_alto && RNA2(r+1,c+1) > thres_bajo
            link(r+1,c+1)=1;
            
            end 
            if RNA2(r+1,c-1) < thres_alto && RNA2(r+1,c-1)> thres_bajo
            link(r+1,c-1)=1;
            
            end

            if RNA2(r-1,c) < thres_alto &&  RNA2(r-1,c) > thres_bajo
            link(r-1,c)=1;
            
            end

            if RNA2(r-1,c+1) < thres_alto && RNA2(r-1,c+1) > thres_bajo 
            link(r-1,c+1)=1;
            
            end

            if RNA2(r-1,c-1) < thres_alto && RNA2(r-1,c-1) > thres_bajo 
            link(r-1,c-1)=1;
            
            end 

            if RNA2(r, c +1) < thres_alto  && RNA2(r, c +1) > thres_bajo 
            link(r, c +1)=1;
            
            end

            if RNA2(r,c-1) < thres_alto &&RNA2(r,c-1) > thres_bajo 
            link(r, c-1)=1;
            
            end 

            if RNA2(r+1,c) < thres_alto &&  RNA2(r+1,c) > thres_bajo
            link(r+1,c)=1;
            
            end 
        end 
    end 
end 
    link_domain2=link;



function out = isboundary(matrix) 
[r,c] = find(matrix);
out=zeros(size(matrix));

for h=1:length(r)
    i=r(h);
    j=c(h);
    if matrix(i,j+1) ==0 || matrix(i,j-1)==0 || matrix (i-1,j)==0 ||matrix(i+1, j)==0
        out(i,j)=1;
    end
end

            
function [RNA,FileName] = func_load(data_directory,AllFileNames,fileindex,resize_factor)

    FileName=AllFileNames(fileindex).name;
    framedir=strcat(data_directory,'\',FileName);   
    RNA=read_RemoveHeader(framedir);
    RNA= imresize(RNA,[resize_factor,resize_factor],'bicubic');

function [RNA,FileName] = func_load_jpg(data_directory,AllFileNames,fileindex,resize_factor)

    FileName=AllFileNames(fileindex).name;
    framedir=strcat(data_directory,'\',FileName);   
    RNA=imread(framedir);
    RNA= imresize(RNA,[resize_factor,resize_factor],'bicubic');





% --- Executes on button press in reset_colors.
function reset_colors_Callback(hObject, eventdata, handles)
% hObject    handle to reset_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.RNAjpg) == 0 
    axes(handles.axes4)
    cla reset
    imshow(handles.RNAjpg)
    
else
    axes(handles.axes4)
    cla reset
    imshow(handles.RNA)
end 

guidata(hObject, handles);



