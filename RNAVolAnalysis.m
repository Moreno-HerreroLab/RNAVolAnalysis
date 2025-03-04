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
% Last Modified by GUIDE v2.5 02-Mar-2025 17:31:35
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

data_directory = 'C:\Users\evatr\OneDrive - UAM\Documentos\MorenoHerreroLab\Projects_tesis\CONCR\CONCR data sorted\128 dataset\txt';
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
handles.color='b';
handles.Remove_Smaller_Than=150;

handles.resize_factor=512;

set(handles.pixel_size,'String', num2str(handles.resize_factor));


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


% --- Executes on button press in browser. Browses a data_path for the
% images folder.

function browser_Callback(~, ~, handles)

input_pathname=uigetdir();
set(handles.data_path,'String', num2str(input_pathname));



function data_path_Callback(hObject, ~, ~)

chckfldr=get(hObject,'String');
chckfldr= exist(chckfldr,'dir');

if chckfldr ~= 7,  msgbox('The provided directory is not valid. Please browse or insert a valid one','ERROR', 'error')
     return; end


% --- Executes on button press in load_first. Load the first image file in
% the folder on data_path 
function load_first_Callback(hObject, ~, handles)

h = warndlg('You are loading a new folder. Please check the number of nucleotides of this sample is correct! Otherwise it will return incorrect Cumulative Sums', 'Warning');
uiwait(h);

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
    [handles.RNA, handles.FileName,handles.factor]= func_load(handles.data_directory,handles.AllFileNames, handles.fileindex,handles.resize_factor);
    isempty(handles.JPGFileNames);
    handles.RNAjpg=[];
    if isempty(handles.JPGFileNames) == 0
        [handles.RNAjpg, handles.JPGFileName]= func_load_jpg(handles.data_directory,handles.JPGFileNames, handles.fileindex,handles.resize_factor);
        handles.RNAjpg_tosave=handles.RNAjpg;

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


% --- Executes on button press in load_next. After loading the first image
% (load_first) it load the next images. 
function load_next_Callback(hObject, ~, handles)

handles.order=0;
if handles.fileindex == handles.chckflder
    msgbox('No next frame','ERROR', 'error')
else
    handles.fileindex=handles.fileindex+1;
    [handles.RNA, handles.FileName,handles.factor]= func_load(handles.data_directory,handles.AllFileNames, handles.fileindex,handles.resize_factor);
    handles.RNAjpg=[];
    set(handles.frame_name,'String', num2str(handles.FileName));
    if isempty(handles.JPGFileNames) == 0
        [handles.RNAjpg, handles.JPGFileName]= func_load_jpg(handles.data_directory,handles.JPGFileNames, handles.fileindex,handles.resize_factor);
        handles.RNAjpg_tosave=handles.RNAjpg;
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
    [handles.RNA, handles.FileName,handles.factor]= func_load(handles.data_directory,handles.AllFileNames, handles.fileindex,handles.resize_factor);
    handles.RNAjpg=[];
    set(handles.frame_name,'String', num2str(handles.FileName));
    if isempty(handles.JPGFileNames) == false
        [handles.RNAjpg, handles.JPGFileName]= func_load_jpg(handles.data_directory,handles.JPGFileNames, handles.fileindex,handles.resize_factor);
        handles.RNAjpg_tosave=handles.RNAjpg;
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



% --- Executes on button press in smooth_button. Applies a gaussian smooth
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



% --- Executes on slider movement. When moved it chages the threshold of
% the molecule 
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


%set(handles.uitable1,'Data',handles.Volume')
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
handles.RNA_molecule=apply_threshold(handles);

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

handles.data = str2double(s);
handles.thres=[handles.thres_array(1) handles.thres_array(end)];
handles.infor = string(handles.info_sort_array);

axes(handles.axes4);

hold on
aux=zeros(size(handles.RNA));
aux(handles.labelMatrix == handles.row)=1;
border=bwperim(aux);
[xx,yy]=find(border);
plot(yy,xx,'.k','MarkerSize',10);
% 
% %s=scatter(y,x,'filled','SizeData',1);
% s=scatter(y,x,'filled','SizeData',1,'MarkerFaceColor',handles.color);
% s.MarkerFaceAlpha = .3;
se = strel('disk',1);
aux=zeros(size(handles.RNAjpg));
aux(border)=1;
aux_dilated = imdilate(aux,se);
[xb,yb]=find(aux_dilated);
for i=1:length(xb)
handles.RNAjpg_tosave(xb(i),yb(i),1)=0;
handles.RNAjpg_tosave(xb(i),yb(i),2)=0;
handles.RNAjpg_tosave(xb(i),yb(i),3)=0;
end 

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

function guardarLinea(filename, celda)
    fid = fopen(filename, 'a'); % 'a' para agregar sin sobrescribir
    if fid == -1
        error('No se pudo abrir el archivo.');
    end

    j = 1;
    while j<length(celda)
        fprintf(fid, '%s,', celda{j});
        datos = celda{j+1};
        if isnumeric(datos)  % Si los datos son numéricos, conviértelos a texto
            fprintf(fid, '%g,', datos(1:end)); % Escribir los valores separados por comas
            
        elseif ischar(datos) || isstring(datos) % Si es texto
            fprintf(fid, '%s,', datos);
        else
            error('Formato de datos no soportado.');
        end
        
        j= j+2;
    end 
    fprintf(fid, '\n'); % Agregar salto de línea al final

    fclose(fid);



% --- Executes on button press in save_molecule.
function save_molecule_Callback(hObject, ~, handles)
% hObject    handle to save_molecule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

st=handles.FileName;
ix=find(st=='.');
%num=str2double(st(1:ix-1));



%First we calculate wich nucleotides are in a domain:
%------------------------------------------------------------------------

N_nucleotides=handles.N_nucleotides;

nucleotides_arr=linspace(1,N_nucleotides,N_nucleotides); %array with nucleotides
total_vol=sum(handles.vol_sort_array);
normalized_volumes=handles.vol_sort_array./total_vol;
nucleotides=normalized_volumes.*N_nucleotides;
cumu_nucleotides=round(cumsum(nucleotides),0);
probabilities=zeros(1,N_nucleotides);
handles.infor
if handles.infor(1)=="Domain" 
        probabilities(nucleotides_arr<=cumu_nucleotides(1))=1;
end 
for i=2:length(cumu_nucleotides)
    if handles.infor(i)=="Domain" 
        probabilities(nucleotides_arr<=cumu_nucleotides(i) & nucleotides_arr>cumu_nucleotides(i-1))=1;
    end 
end 

%------------------------------------------------------------------------

total_cell = {st, handles.vol_sort_array,'info', handles.infor,'thresholds',handles.thres_array,'background noise', handles.Noise};

guardarLinea('CumuSums.csv',{st cumu_nucleotides})
guardarLinea('Results.csv', total_cell)
guardarLinea('Probabilities.csv',{st probabilities})

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



function pixel_size_Callback(hObject, ~, handles)
% hObject    handle to pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixel_size as text
%        str2double(get(hObject,'String')) returns contents of pixel_size as a double

handles.resize_factor=str2double(get(handles.pixel_size,'String'));
%handles.factor=handles.size_nm/handles.resize_factor;


guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function pixel_size_CreateFcn(hObject, ~, ~)
% hObject    handle to pixel_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function [RNA1, xAmplitude]=read_RemoveHeader(path)
fid1 = fopen(path,'r');
linedata = '0';
lineNumber = 0;
mat=[];
xAmplitude = [];

while (ischar(linedata))
    linedata = fgetl(fid1);
    lineNumber = lineNumber +1;
    if lineNumber <= 5 & contains(linedata, 'X Amplitude:')
        xAmplitude = sscanf(linedata, 'X Amplitude: %f nm');
    end
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

if get(handles.checkbox2,'Value')==1
    RNA3=handles.ROI;
else
    RNA3=handles.RNA2;
end 

mask=RNA3>handles.threshold;
RNA3(mask==0)=0;


axes(handles.axes2)
cla reset
imshow(RNA3);


function [linkers,objects]=main_run(handles)

handles.thres_array=sort(handles.thres_array);

n_steps=length(handles.thres_array);
linker_map=zeros(size(handles.RNA2));
object_map=zeros(size(handles.RNA2));

if get(handles.checkbox2,'Value')==1
    handles.RNA2=handles.ROI;
end 



for i=1:n_steps-1
    low=handles.thres_array(i);
    high=handles.thres_array(i+1);
    
    %guardamos la lista de pixeles de cada isla con el thres alto
    maskhigh=zeros(size(handles.RNA2));
    maskhigh(handles.RNA2 > high) = 1; 
    %maskhigh=imclearborder(maskhigh);
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
    %de pixeles de las regiones que van a unirse (thres alto) 
    
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
            %guardamos las islas que se han unido
             prev_objects_idx=find(tf==1);
             
            %objetos se han unido en uno solo
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
                    [~,mask]=find_8_conn2(ep_loc(ep_i),handles.RNA2,0);
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


link_border=bwperim(linker_map);
[x, y] = find(link_border);
plot(y,x,'r.','LineWidth',10);

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
if get(handles.checkbox2,'Value')==1
    RNA2=handles.ROI;
else
    RNA2=handles.RNA2;
end 

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
for i=1:total_dom
    [x,y]=find(handles.labelMatrix == i);   
    handles.RNA_fig=insertText(handles.RNA_fig, [mean(y) mean(x)], num2str(i),'FontSize',10);
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



            
function [RNA,FileName,factor] = func_load(data_directory,AllFileNames,fileindex,resize_factor)

    FileName=AllFileNames(fileindex).name;
    framedir=strcat(data_directory,'\',FileName);   
    [RNA,XAmplitude]=read_RemoveHeader(framedir);
    factor = XAmplitude/resize_factor;
    RNA= imresize(RNA,[resize_factor,resize_factor],'bicubic');

function [RNA,FileName] = func_load_jpg(data_directory,AllFileNames,fileindex,resize_factor)

    FileName=AllFileNames(fileindex).name;
    framedir=strcat(data_directory,'\',FileName);   
    RNA=imread(framedir);
    RNA= imresize(RNA,[resize_factor,resize_factor],'bicubic');
    



function reset_colors_Callback(hObject, ~, handles)
% hObject    handle to reset_colors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.RNAjpg) == 0 
    axes(handles.axes4)
    cla reset
    imshow(handles.RNAjpg)
    handles.RNAjpg_tosave=handles.RNAjpg;
else
    axes(handles.axes4)
    cla reset
    imshow(handles.RNA)
end 

guidata(hObject, handles);



function adjust_Callback(hObject, ~, handles)


axes(handles.axes2);
imcontrast

guidata(hObject, handles);



function filter_1_Callback(hObject, ~, handles)


message = sprintf('Left click and hold to begin drawing.\nSimply lift the mouse button to finish');
uiwait(msgbox(message));
hFH = drawfreehand(handles.axes2);

binaryMask = createMask(hFH);
handles.ROI=handles.RNA2;
handles.ROI(binaryMask==0)=0;

axes(handles.axes6);
imshow(handles.ROI)


guidata(hObject, handles);



function Add_Callback(hObject, ~, handles)
% hObject    handle to Add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


vol1=handles.Volume(handles.row).vol;
handles.tot_sum=handles.tot_sum + vol1;
handles.NewReg(handles.labelMatrix==handles.row)=1;

handles.Volume(handles.row).vol=0;
handles.labelMatrix(handles.labelMatrix==handles.row)=0;

guidata(hObject, handles);

% --- Executes on button press in Equal.
function Equal_Callback(hObject, ~, handles)
% hObject    handle to Equal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Volume2(1).vol=handles.tot_sum;
handles.Volume2(1).info='Domain';
handles.labelMatrix2=handles.NewReg;
B=unique(handles.labelMatrix);
B=nonzeros(B);
for i=1:length(B)
    handles.labelMatrix2(handles.labelMatrix==B(i))=i+1;
    handles.Volume2(i+1).vol=handles.Volume(B(i)).vol;
    handles.Volume2(i+1).info=handles.Volume(B(i)).info;
end 

handles.Volume=handles.Volume2;
% set(handles.uitable1,'Data',handles.Volume2')
cell1={handles.Volume2.vol};
cell2={handles.Volume2.info};
handles.TAB=[cell1;cell2];

set(handles.uitable1,'Data',handles.TAB')


axes(handles.axes3);
cla reset
hold on 
handles.RNA_fig=handles.RNA;

for i=1:length(handles.Volume2)
    [x,y]=find(handles.labelMatrix2 == i);   
    handles.RNA_fig=insertText(handles.RNA_fig, [mean(y) mean(x)], num2str(i),'FontSize',10);
end 
imshow(handles.RNA_fig)

for i=1:length(handles.Volume2)
    positt=handles.labelMatrix2 == i;
    aux=zeros(size(handles.RNA));
    aux(positt)=1;
    border=bwperim(aux);
    [xx,yy]=find(border);
    plot(yy,xx,'.r','MarkerSize',3);
    
end 

handles.labelMatrix=handles.labelMatrix2;


handles.Volume2=[];
handles.tot_sum=0;
handles.NewReg=zeros(size(handles.RNA2));

guidata(hObject, handles);


% --- Executes on button press in inizialize_sum.
function inizialize_sum_Callback(hObject, ~, handles)
% hObject    handle to inizialize_sum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Volume2=[];
handles.tot_sum=0;
handles.NewReg=zeros(size(handles.RNA2));

guidata(hObject, handles);





% --- Executes on button press in skeletonize_button.
function skeletonize_button_Callback(hObject, ~, handles)
% hObject    handle to skeletonize_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[m,n]=size(handles.RNA2);
handles.skeleton_image=ones(m,n,3);

handles.RNA_molecule=imclearborder(handles.RNA_molecule);

handles.skeleton = bwmorph(handles.RNA_molecule,'thin', Inf);
handles.skeleton = bwareaopen( handles.skeleton , handles.Remove_Smaller_Than );


se = strel('disk',5);
handles.skeleton_dilated = imdilate(handles.skeleton,se);
[x,y]=find(handles.skeleton_dilated);
for i=1:length(x)
    handles.skeleton_image(x(i),y(i),1)=0;
    handles.skeleton_image(x(i),y(i),2)=0;
    handles.skeleton_image(x(i),y(i),3)=0;
end 
axes(handles.axes6);
cla reset
imshow(handles.skeleton_image);

if isempty(handles.RNAjpg) == 0 
    axes(handles.axes8)
    cla reset
    imshow(handles.RNAjpg)
    
else
    axes(handles.axes8)
    cla reset
    imshow(handles.RNA)
end 
axes(handles.axes8);
hold on
[x,y]=find(handles.skeleton_image == 0);
scatter(y,x,'filled','SizeData',1,'MarkerFaceColor','k');
handles.RNAskel_tosave=handles.skeleton_image;

guidata(hObject, handles);


ep = bwmorph(handles.skeleton,'endpoints');
%Positions of the end-points
[r,c] = find(ep); %Position (r:row, c:colum)
ends_loc = find(ep); %Position (index)
endpoints=[r,c,ends_loc];

bp = bwmorph(handles.skeleton,'branchpoints');
[r1,c1] = find(bp); 
B_loc = find(bp);
branchpoints=[r1,c1,B_loc];

%Inizialite max_length variable
max_length=0;
dd=size(endpoints);

%Here we find the main chain distance 
for k = 1:dd(1)
    D = bwdistgeodesic(handles.skeleton,endpoints(k,2),endpoints(k,1),'quasi-euclidean'); 
    D(isinf(D)) = nan;
    endtoend(k) = max(D(endpoints(:,3))); 
    ind = D == endtoend(k); %Keeps the position of the endpoint to which [c(k),r(k)] distance is max
    [r2,c2]= find(ind); %Keeps the same position in r:row and c:colum 
    index=find(ind);
    if endtoend(k) > max_length && (endtoend(k) ~= inf)
        max_length = endtoend(k);
        ep1(1) = endpoints(k,1); %one endpoint of the longest chain
        ep1(2) = endpoints(k,2);
        ep1(3) = endpoints(k,3);
        ep2(1) = r2; %The second endpoint of the main chain 
        ep2(2) = c2;
        ep2(3) = index;
    end
end

D1 = bwdistgeodesic(handles.skeleton, ep1(2), ep1(1), 'quasi-euclidean'); %endpoint 1
D2 = bwdistgeodesic(handles.skeleton, ep2(2), ep2(1), 'quasi-euclidean'); %endpoint 2
D = D1 + D2;
D = round(D * 8) / 8;
D(isnan(D)) = inf;

mc_lenght=max_length.*handles.factor;


%Now we represent in the image the branches lenght (in red)

if isempty(B_loc) == true
    handles.mc= mc_lenght;
    handles.skeleton_lenght = mc_lenght;
    handles.dob = 1;
    s=regionprops(handles.skeleton,'centroid');
    centroids=s.Centroid;
    axes(handles.axes6);
    hold on;
    plot(centroids(1),centroids(2),'b*')
    text(centroids(1),centroids(2), sprintf('%.2f', mc_lenght), ...
        'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
else
    for k = 1:numel(r)
        matrix_resta=false(size(handles.skeleton));
        D = bwdistgeodesic(handles.skeleton,endpoints(k,2),endpoints(k,1),'quasi-euclidean');
        distanceToBranchPt = min(D(B_loc)); %Distances from the endpoint [c(k),r(k)] to all bp. 
                                        %the shortest parth is the length of the branch.
                                        %Takes into account the starting and ending branch of the main chain      
        branches_length(k) = distanceToBranchPt.*handles.factor;
        matrix_resta(D < distanceToBranchPt)=true;
        s=regionprops(matrix_resta,'centroid');
        centroids=s.Centroid;
        axes(handles.axes6);
        hold on;
        plot(centroids(1),centroids(2),'b*')
        text(centroids(1),centroids(2), sprintf('%.2f', branches_length(k)), ...
            'Color', 'red', 'FontSize', 12, 'FontWeight', 'bold');
    end 

hold off;

%Now we calculate the lenght of the skeleton

skel_length=0;
skeleton_aux=handles.skeleton;
while length(B_loc) >= 2
    branches_length=[];
    matrix_resta=false(size(handles.skeleton));
    for k = 1:numel(r)
            D = bwdistgeodesic(skeleton_aux,endpoints(k,2),endpoints(k,1),'quasi-euclidean');
            distanceToBranchPt = min(D(B_loc)); %Distances from the endpoint [c(k),r(k)] to all bp. 
                                            % the shortest parth is the length of the branch.
                                            % Takes into account the starting and ending branch of the main chain                  
            branches_length(k) = distanceToBranchPt.*handles.factor;
            matrix_resta(D < distanceToBranchPt)=true; %Matrix_resta contains all the branches
    end 
    skeleton_aux=skeleton_aux - matrix_resta; %Skeleton_aux contains a new skeleton (original - branches)
    skeleton_aux=bwmorph(skeleton_aux,'thin');

    %Calculate new ep and bp of the new skeleton 
    ep = bwmorph(skeleton_aux,'endpoints');
    [r,c] = find(ep); %Position (r:row, c:colum)
    ends_loc = find(ep); %Position (index)
    endpoints=[r,c,ends_loc];

    bp = bwmorph(skeleton_aux,'branchpoints');
    %[~,~]=find(bp);
    B_loc = find(bp); %Update the value of B_loc (location of bp) - When number of bp is 1 or 0, the while loop ends
    skel_length = skel_length + sum(branches_length);
end 

core_molecule = skeleton_aux;

%2 cases: 1 bp remains, or no bp remains.
if length(B_loc)==1 %If one bp remains, we erase branches until 0 bp remains:
    while length(B_loc)==1
        D = bwdistgeodesic(skeleton_aux,endpoints(1,2),endpoints(1,1),'quasi-euclidean');
        distanceToBranchPt = min(D(B_loc)); %Distances from the endpoint [c(k),r(k)] to all bp. 
                                        % the shortest parth is the length of the branch.
                                        % Takes into account the starting and ending branch of the main chain                  
        matrix_resta(D < distanceToBranchPt)=true; %Matrix_resta contains one branches
        skeleton_aux=skeleton_aux - matrix_resta; 
        skeleton_aux=bwmorph(skeleton_aux,'thin');
        skel_length = skel_length + distanceToBranchPt.*handles.factor;
        %Calculate new ep and bp of the new skeleton 
        ep = bwmorph(skeleton_aux,'endpoints');
        [r,c] = find(ep); %Position (r:row, c:colum)
        ends_loc = find(ep); %Position (index)
        endpoints=[r,c,ends_loc];
    
        bp = bwmorph(skeleton_aux,'branchpoints');
        B_loc = find(bp); %Update the value of B_loc (location of bp) - When number of bp is 1 or 0, the while loop ends
    end 
else %If no bp remains, we add the last segment:
    D = bwdistgeodesic(skeleton_aux,endpoints(1,2),endpoints(1,1),'quasi-euclidean');
    D(D==0)=inf;
    distanceToBranchPt=min(D(ends_loc));
    matrix_resta(D < distanceToBranchPt)=true;
    if distanceToBranchPt ~= inf
        skel_length = skel_length + distanceToBranchPt.*handles.factor;
    end 
end 

%Final values of skeleton lenghts:
handles.mc= mc_lenght;
handles.skeleton_lenght = skel_length;
handles.dob = handles.skeleton_lenght/mc_lenght;

% We represent in the image the remaining values (segment between
% branchpoints mainly)

bp = bwmorph(handles.skeleton,'branchpoints');
bp_locs=find(bp);

[~,mask]=find_8_conn2(bp_locs,handles.skeleton,1);

%Now the skeleton is divided in segments 
I_segmented = core_molecule & ~mask;
[L, num] = bwlabel(I_segmented);

for k = 1:num
    matrix_resta=false(size(handles.skeleton));
    branch_mask = L == k;
    ep = bwmorph(branch_mask,'endpoints');
    [r,c]=find(ep);
    ends_loc = find(ep);
    %Find the distance from every point in the branch to the nearest endpoint
    D = bwdistgeodesic(branch_mask,c(1) ,r(1), 'quasi-euclidean');
    distance = max(D(ends_loc));  % Max distance will be the length of the branch

    matrix_resta(D < distance )=true;
    s=regionprops(matrix_resta,'centroid');
    centroids=s.Centroid;
    axes(handles.axes6);
    hold on;
    plot(centroids(1),centroids(2),'b*')
    text(centroids(1),centroids(2), sprintf('%.2f', distance.*handles.factor), ...
        'Color', 'blue', 'FontSize', 12, 'FontWeight', 'bold');
end
hold off
end 

etiquetas = bwlabel(handles.RNA_molecule);
props = regionprops(etiquetas, 'Area');
[~, indice_max_area] = max([props.Area]);
RNA3_main = (etiquetas == indice_max_area);

props1= regionprops(RNA3_main, 'MajorAxisLength',"MaxFeretProperties", "MinFeretProperties");
props = regionprops(RNA3_main, 'BoundingBox');

% Calcular el diámetro mínimo basado en el cuadro delimitador mínimo
ancho = props.BoundingBox(3);
alto = props.BoundingBox(4);
handles.diametro_minimo = sqrt(ancho^2 + alto^2)*handles.factor;

handles.MajorAxisLength =props1.MajorAxisLength*handles.factor;
handles.MaxFeret= props1.MaxFeretDiameter*handles.factor;
handles.MinFeret= props1.MinFeretDiameter*handles.factor;


guidata(hObject, handles);




% --- Executes on button press in saveSkel.
function saveSkel_Callback(hObject, ~, handles)
% hObject    handle to saveSkel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


st=handles.FileName;


total_cell ={st,handles.mol_thres,'Skeleton Length',handles.skeleton_lenght,'Main Chain Length',handles.mc,'MajorAxis Lenght',handles.MajorAxisLength,'Smallest diameter',handles.diametro_minimo,'MinFeret',handles.MinFeret,'MaxFeret',handles.MaxFeret};
guardarLinea('Results_Skels.csv',total_cell)
guidata(hObject, handles);

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
% if get(hObject,'Value') ==0
%     handles.ROI=0;
% end 
guidata(hObject, handles);


% --- Executes on button press in red_butt.
function red_butt_Callback(hObject, ~, handles)
% hObject    handle to red_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes6);  % Set the current axes
[x, y] = ginput(1);   % Get the user's click
hold on;
plot(x, y, 'ro', 'MarkerSize', 2,'LineWidth',5);  % Plot a blue dot
hold off;


se=strel('disk',8);
aux=zeros(size(handles.RNAskel_tosave));
aux(round(y,0),round(x,0))=1;
aux_dilated = imdilate(aux,se);
[xb,yb]=find(aux_dilated);
for i=1:length(xb)
handles.RNAskel_tosave(xb(i),yb(i),1)=1;
handles.RNAskel_tosave(xb(i),yb(i),2)=0;
handles.RNAskel_tosave(xb(i),yb(i),3)=0;
end 


guidata(hObject, handles);

% --- Executes on button press in blue_butt.
function blue_butt_Callback(hObject, ~, handles)
% hObject    handle to blue_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes6);  % Set the current axes
[x, y] = ginput(1);   % Get the user's click
hold on;
plot(x, y, 'bo', 'MarkerSize', 2,'LineWidth',5);  % Plot a blue dot
hold off;

se=strel('disk',8);
aux=zeros(size(handles.RNAskel_tosave));
aux(round(y,0),round(x,0))=1;
aux_dilated = imdilate(aux,se);
[xb,yb]=find(aux_dilated);

for i=1:length(xb)
handles.RNAskel_tosave(xb(i),yb(i),1)=0;
handles.RNAskel_tosave(xb(i),yb(i),2)=0;
handles.RNAskel_tosave(xb(i),yb(i),3)=1;
end 



guidata(hObject, handles);






function [extended_bp_locations,mask]=find_8_conn2(b_locs,RNA, keep_bp)

mask=zeros(size(RNA));
mask(b_locs)=1;

[r,c]=find(mask);

if keep_bp == 0
    mask=zeros(size(RNA));
end 

for i=1:numel(r)
    mask(r(i)+1,c(i))=1;
    mask(r(i)-1,c(i))=1;
    mask(r(i),c(i)+1)=1;
    mask(r(i),c(i)-1)=1;
    mask(r(i)-1,c(i)-1)=1;
    mask(r(i)+1,c(i)+1)=1;
    mask(r(i)+1,c(i)-1)=1;
    mask(r(i)-1,c(i)+1)=1;
end 

    
mask=imclearborder(mask);
extended_bp_locations=find(mask);


function [link_domain2]=extend_border_link(link,thres_alto,thres_bajo,RNA2)

[r_total,c_total]=size(RNA2);

for h=1:3
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



% --- Executes on button press in browse_threshods.
function browse_threshods_Callback(hObject, ~, handles)
% hObject    handle to browse_threshods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[filename, filepath] = uigetfile('*.txt', 'Selecciona un archivo');

% Combinar el nombre del archivo y la ubicación para obtener la ruta completa
nombre_archivo = fullfile(filepath, filename);
% Abrir el archivo para lectura
fid = fopen(nombre_archivo, 'r');

% Inicializar el vector para almacenar los valores de la tercera columna
thresholds = [];
file_ids=[];
% Leer el archivo línea por línea
while ~feof(fid)
    % Leer una línea del archivo
    linea = fgetl(fid);
    
    % Dividir la línea en partes utilizando ';' como delimitador
    partes = strsplit(linea, ';');
    
    % Obtener el valor de la tercera columna y convertirlo a número
    valor = str2double(partes{3});

    file_id = str2double(partes{1});
    
    % Agregar el valor al vector thresholds
    thresholds = [thresholds; valor];
    file_ids=[file_ids;file_id];
end

% Cerrar el archivo
fclose(fid);

handles.All_Thresholds=thresholds;
handles.Files_Ids=file_ids;


guidata(hObject, handles);


% --- Executes on button press in Automatic_Load.
function Automatic_Load_Callback(hObject, ~, handles)
% hObject    handle to Automatic_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.order=0;
handles.data_directory=get(handles.data_path,'String');
CurrentFolder=pwd;
handles.fileindex=1;
handles.fileindex_load=1;
cd(handles.data_directory);

handles.AllFileNames=dir('*.txt');
cd(CurrentFolder)
handles.chckflder=length(handles.AllFileNames);

if handles.chckflder == 0
    msgbox('The provided input folder is empty.','ERROR', 'error')
else 
    MajorAxisLengths= [];
    diametros_minimos =[];
    MaxFeret_all=[];
    MinFeret_all=[];
    for h=1:length(handles.All_Thresholds)
    [handles.RNA, handles.FileName, handles.factor]= func_load(handles.data_directory,handles.AllFileNames, handles.fileindex_load,handles.resize_factor);
    st=handles.FileName;
    ix=find(st=='.');
    num=str2double(st(1:ix-1));
    while num ~= handles.Files_Ids(handles.fileindex)
        handles.fileindex_load= handles.fileindex_load+1;
        [handles.RNA, handles.FileName,handles.factor]= func_load(handles.data_directory,handles.AllFileNames, handles.fileindex_load,handles.resize_factor);
        st=handles.FileName;
        ix=find(st=='.');
        num=str2double(st(1:ix-1));
    end 
    [MajorAxisLength,diametro_minimo,MaxFeret,MinFeret,skeleton_matrix]=automatic_run(handles.RNA,handles.fileindex,handles.All_Thresholds,handles.axes1); 
    
    axes(handles.axes1)
    imshow(skeleton_matrix)
    pause
    name=[st(1:ix-1),'.txt'];
    writematrix(skeleton_matrix,  name,'Delimiter','tab');
    MajorAxisLengths = [MajorAxisLengths;MajorAxisLength];
    diametros_minimos =[diametros_minimos; diametro_minimo];
    MaxFeret_all=[MaxFeret_all;MaxFeret];
    MinFeret_all=[MinFeret_all;MinFeret];
    %pause
    handles.fileindex = handles.fileindex+1;
    handles.fileindex_load = handles.fileindex_load+1;
    handles.order=0;
    end 
end 
% 
% cell2={MajorAxisLengths};
% cell3={diametros_minimos};
% 
% handles.Results_Props=[cell2;cell3];
% writecell(handles.Results_Props,'EnclosingCircleData.dat');

fid = fopen('Results_Props.txt', 'w');
fprintf(fid, 'MajorAxisLengths\tdiametro_minimo\tMaxFeret\tMinFeret\n'); % Encabezado
fprintf(fid, '%f\t%f\t%f\t%f\n', [MajorAxisLengths, diametros_minimos,MaxFeret_all,MinFeret_all].'); % Datos
fclose(fid);

guidata(hObject, handles);


function [MajorAxisLength,diametro_minimo,MaxFeret,MinFeret,skeleton_matrix] = automatic_run (RNA,fileindex,All_Thresholds,axess)

handles.RNA2=imgaussfilt(RNA,5);
RNA3=handles.RNA2;
mask=RNA3>All_Thresholds(fileindex);
RNA3(mask==0)=0;

etiquetas = bwlabel(RNA3);

% Calcular el área de cada región etiquetada
props = regionprops(etiquetas, 'Area');
% Encontrar la etiqueta de la región con el área más grande
[~, indice_max_area] = max([props.Area]);

% Crear una nueva imagen binaria con solo la región más grande
RNA3_main = (etiquetas == indice_max_area);
% 
% axes(axess);
% cla reset
% imshow(RNA3_main);

props1= regionprops(RNA3_main, 'MajorAxisLength',"MaxFeretProperties", "MinFeretProperties");
props = regionprops(RNA3_main, 'BoundingBox');

% Calcular el diámetro mínimo basado en el cuadro delimitador mínimo
ancho = props.BoundingBox(3);
alto = props.BoundingBox(4);
diametro_minimo = sqrt(ancho^2 + alto^2);

MajorAxisLength =props1.MajorAxisLength;
MaxFeret= props1.MaxFeretDiameter;
MinFeret= props1.MinFeretDiameter;


% ahora calcular skeleton e importar matix con el numero.txt

skeleton_matrix = bwmorph(RNA3_main,'thin', Inf);




function edit12_Callback(hObject, ~, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


handles.Remove_Smaller_Than= str2double(get(hObject,'String')) ;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, ~, ~)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N_nucleotides_Callback(hObject, eventdata, handles)
% hObject    handle to N_nucleotides (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_nucleotides as text
%        str2double(get(hObject,'String')) returns contents of N_nucleotides as a double

handles.N_nucleotides= str2double(get(hObject,'String')) ;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function N_nucleotides_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_nucleotides (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


