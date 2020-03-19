% matlab gui to annotate stroke phantom
% author: michael manhart | michael.manhart@informatik.uni-erlangen.de
%         pattern recognition lab, university of erlangen-nuremberg
% last change: 16.01.2012

function varargout = strokecreator(varargin)
% STROKECREATOR MATLAB code for strokecreator.fig
%      STROKECREATOR, by itself, creates a new STROKECREATOR or raises the existing
%      singleton*.
%
%      H = STROKECREATOR returns the handle to a new STROKECREATOR or the handle to
%      the existing singleton*.
%
%      STROKECREATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STROKECREATOR.M with the given input arguments.
%
%      STROKECREATOR('Property','Value',...) creates a new STROKECREATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before strokecreator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to strokecreator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help strokecreator

% Last Modified by GUIDE v2.5 19-Jan-2012 09:50:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @strokecreator_OpeningFcn, ...
                   'gui_OutputFcn',  @strokecreator_OutputFcn, ...
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


% --- Executes just before strokecreator is made visible.
function strokecreator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to strokecreator (see VARARGIN)

% Choose default command line output for strokecreator
handles.output = hObject;

% initialization
[handles.brain_o handles.brain handles.mrbrain] = loadBrain();
handles.sliceNum = round(size(handles.brain,3)/2);
handles.axis = [1 size(handles.brain,1) 1 size(handles.brain,2) 1 size(handles.brain,3)];
handles.showMR = 0;
handles.annotationMode = 1;     % annotationMode: 1 = freehand, 2 = 2D Ball, 3 = 3D Ball
handles.viewDim = 'z';
set(handles.slider_slice, 'value', handles.sliceNum, 'min', handles.axis(5), 'max', handles.axis(6), 'sliderstep', [1/(handles.axis(6)- handles.axis(5)) 1/(handles.axis(6)- handles.axis(5))]);
set(handles.text_slice, 'String', num2str(handles.sliceNum));  
set(handles.uipanel_annotationMode,'SelectionChangeFcn',@uipanel_annotationMode_SelectionChangeFcn);
showBrain(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes strokecreator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = strokecreator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_slice_Callback(hObject, eventdata, handles)
% hObject    handle to slider_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sliceNum = round(get(hObject,'Value'));
showBrain(handles);
guidata(hObject,handles);
set(handles.text_slice, 'String', num2str(handles.sliceNum));  
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_slice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_slice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popup_selectDimension.
function popup_selectDimension_Callback(hObject, eventdata, handles)
% hObject    handle to popup_selectDimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
handles.viewDim = contents{get(hObject,'Value')};
switch handles.viewDim
    case 'z'
        handles.sliceNum = round(size(handles.brain,3)/2);       
    case 'y'
        handles.sliceNum = round(size(handles.brain,2)/2);
    case 'x'
        handles.sliceNum = round(size(handles.brain,3)/2);        
end
set(handles.text_slice, 'String', num2str(handles.sliceNum));  
switch handles.viewDim
    case 'z'
        set(handles.slider_slice, 'value', handles.sliceNum, 'min', handles.axis(5), 'max', handles.axis(6), 'sliderstep', [1/(handles.axis(6)- handles.axis(5)) 1/(handles.axis(6)- handles.axis(5))])
    case 'y'
        set(handles.slider_slice, 'value', handles.sliceNum, 'min', handles.axis(3), 'max', handles.axis(4), 'sliderstep', [1/(handles.axis(4)- handles.axis(3)) 1/(handles.axis(4)- handles.axis(3))])
    case 'x'
        set(handles.slider_slice, 'value', handles.sliceNum, 'min', handles.axis(1), 'max', handles.axis(2), 'sliderstep', [1/(handles.axis(2)- handles.axis(1)) 1/(handles.axis(2)- handles.axis(1))])
end
guidata(hObject,handles);
showBrain(handles);
% Hints: contents = cellstr(get(hObject,'String')) returns popup_selectDimension contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_selectDimension


% --- Executes during object creation, after setting all properties.
function popup_selectDimension_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_selectDimension (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_reducedPerfusion.
function button_reducedPerfusion_Callback(hObject, eventdata, handles)
% hObject    handle to button_reducedPerfusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = annotateStroke(handles, 0);
guidata(hObject,handles);
showBrain(handles);

% --- Executes on button press in button_severlyReducedPerfusion.
function button_severlyReducedPerfusion_Callback(hObject, eventdata, handles)
% hObject    handle to button_severlyReducedPerfusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = annotateStroke(handles, 1);
guidata(hObject,handles);
showBrain(handles);

% --- Executes on button press in pushbutton_healthyeval.
function pushbutton_healthyeval_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_healthyeval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = annotateStroke(handles, 2);
guidata(hObject,handles);
showBrain(handles);

% --- Executes on button press in button_clearStroke.
function button_clearStroke_Callback(hObject, eventdata, handles)
% hObject    handle to button_clearStroke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.viewDim ~= 'z'
    msgbox('Annotation only in z viewing dimension possible!','Error');
    return;
end
h = imfreehand(handles.axis_brain);
mask = h.createMask;

% classes of tissue
GM = 1; WM = 2; GMR = 3; GMSR = 4; WMR = 5; WMSR = 6; AI = 7; VO = 8; CSF = 9; HEGM = 10; HEWM = 11;

brainslice = squeeze(handles.brain(:,:,handles.sliceNum));
brainslice_o = squeeze(handles.brain_o(:,:,handles.sliceNum));

brainslice( (brainslice_o == GM) & mask) = GM;
brainslice( (brainslice_o == WM) & mask) = WM;
handles.brain(:,:,handles.sliceNum) = brainslice;

delete(h);
guidata(hObject,handles);
showBrain(handles);

% --- Executes on button press in button_Save.
function button_Save_Callback(hObject, eventdata, handles)
% hObject    handle to button_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
brain = handles.brain;
save('brain.mat','brain');

% --- Executes on button press in button_copyLastSlide.
function button_copyLastSlide_Callback(hObject, eventdata, handles)
% hObject    handle to button_copyLastSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% classes of tissue
if handles.viewDim ~= 'z'
    msgbox('Annotation only in z viewing dimension possible!','Error');
    return;
end
GM = 1; WM = 2; GMR = 3; GMSR = 4; WMR = 5; WMSR = 6; AI = 7; VO = 8; CSF = 9; HEGM = 10; HEWM = 11;
lastSlice = handles.sliceNum-1;
if(lastSlice < 1)
    return;
end

brainslice = squeeze(handles.brain(:,:,handles.sliceNum));
brainslice_o = squeeze(handles.brain_o(:,:,handles.sliceNum));
brainslice_last = squeeze(handles.brain(:,:,lastSlice));

annotated = (brainslice_last == HEGM) | (brainslice_last == HEWM);
brainslice( (brainslice_o == GM) & annotated ) = HEGM;
brainslice( (brainslice_o == WM) & annotated ) = HEWM;
annotated = (brainslice_last == GMR) | (brainslice_last == WMR);
brainslice( (brainslice_o == GM) & annotated ) = GMR;
brainslice( (brainslice_o == WM) & annotated ) = WMR;
annotated = (brainslice_last == GMSR) | (brainslice_last == WMSR);
brainslice( (brainslice_o == GM) & annotated ) = GMSR;
brainslice( (brainslice_o == WM) & annotated ) = WMSR;
handles.brain(:,:,handles.sliceNum) = brainslice;

guidata(hObject,handles);
showBrain(handles);

% --- Executes on button press in button_copyNextSlide.
function button_copyNextSlide_Callback(hObject, eventdata, handles)
% hObject    handle to button_copyNextSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% classes of tissue
if handles.viewDim ~= 'z'
    msgbox('Annotation only in z viewing dimension possible!','Error');
    return;
end
GM = 1; WM = 2; GMR = 3; GMSR = 4; WMR = 5; WMSR = 6; AI = 7; VO = 8; CSF = 9; HEGM = 10; HEWM = 11;
nextSlice = handles.sliceNum+1;
maxSlice = size(handles.brain,3);
if(nextSlice > maxSlice)
    return;
end

brainslice = squeeze(handles.brain(:,:,handles.sliceNum));
brainslice_o = squeeze(handles.brain(:,:,handles.sliceNum));
brainslice_next = squeeze(handles.brain(:,:,nextSlice));

annotated = (brainslice_next == HEGM) | (brainslice_next == HEWM);
brainslice( (brainslice_o == GM) & annotated ) = HEGM;
brainslice( (brainslice_o == WM) & annotated ) = HEWM;
annotated = (brainslice_next == GMR) | (brainslice_next == WMR);
brainslice( (brainslice_o == GM) & annotated ) = GMR;
brainslice( (brainslice_o == WM) & annotated ) = WMR;
annotated = (brainslice_next == GMSR) | (brainslice_next == WMSR);
brainslice( (brainslice_o == GM) & annotated ) = GMSR;
brainslice( (brainslice_o == WM) & annotated ) = WMSR;
handles.brain(:,:,handles.sliceNum) = brainslice;

guidata(hObject,handles);
showBrain(handles);

% --- Executes on button press in button_restoreBrain.
function button_restoreBrain_Callback(hObject, eventdata, handles)
% hObject    handle to button_restoreBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
answer = questdlg('Are you sure?', 'Restore Brain', 'No');
if ~strcmp(answer,'Yes')
    return;
end
handles.brain = handles.brain_o;
guidata(hObject,handles);
showBrain(handles);

% --- Executes on button press in checkbox_showMR.
function checkbox_showMR_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showMR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.showMR = get(hObject,'Value');
guidata(hObject,handles);
showBrain(handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox_showMR

function uipanel_annotationMode_SelectionChangeFcn(hObject, eventdata)
%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 

switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton_freehand'
        handles.annotationMode = 1;
    case 'radiobutton_2DBall'
        handles.annotationMode = 2;
    case 'radiobutton_3DBall'
        handles.annotationMode = 3;
end
%updates the handles structure
guidata(hObject, handles);

function edit_zsize3DEllipse_Callback(hObject, eventdata, handles)
% hObject    handle to edit_zsize3DEllipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_zsize3DEllipse as text
%        str2double(get(hObject,'String')) returns contents of edit_zsize3DEllipse as a double


% --- Executes during object creation, after setting all properties.
function edit_zsize3DEllipse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_zsize3DEllipse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS
% END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS
% END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS % END CALLBACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function showBrain(handles)

if handles.showMR
  switch handles.viewDim
    case 'z'
        brainslice = squeeze(handles.mrbrain(:,:,handles.sliceNum));
    case 'y'
        brainslice = squeeze(handles.mrbrain(:,handles.sliceNum,:));
    case 'x'
        brainslice = squeeze(handles.mrbrain(handles.sliceNum,:,:));
  end 
  iptsetpref('ImshowAxesVisible','off');
  imshow(double(brainslice), [0 490], 'Parent', handles.axis_brain);
  colormap gray;
else
    switch handles.viewDim
        case 'z'
            brainslice = squeeze(handles.brain(:,:,handles.sliceNum));
        case 'y'
            brainslice = squeeze(handles.brain(:,handles.sliceNum,:));
        case 'x'
            brainslice = squeeze(handles.brain(handles.sliceNum,:,:));
    end    
    iptsetpref('ImshowAxesVisible','off');
    imshow(double(brainslice), [0 11], 'Parent', handles.axis_brain);
    colormap(tissueclasses);
end

function handles=annotateStroke(handles, type)
if handles.viewDim ~= 'z'
    msgbox('Annotation only in z viewing dimension possible!','Error');
    return;
end

switch handles.annotationMode
    case 1
        h = imfreehand(handles.axis_brain);
    case {2,3}
        h = imellipse(handles.axis_brain);
        position = h.getPosition;
        xpos = position(1) + position(3)/2;
        ypos = position(2) + position(4)/2;  
        zpos = handles.sliceNum;
        if handles.annotationMode == 2
            zsize = 1;
        else
            zsize = str2double(get(handles.edit_zsize3DEllipse,'String'));
            if isnan(zsize)
                msgbox('Invallid z size of 3D ellipse!','Error');
                return;                
            end       
        end
        disp(['ellipse: xpos: ' num2str(xpos) ' ypos: ' num2str(ypos) ' zpos: ' num2str(zpos)]);
        disp(['        xsize: ' num2str(position(3)) ' mm ysize: ' ...
         num2str(position(4)) ' mm zsize:' num2str(zsize) ' mm' ...
         'volume: ' num2str(position(3)*position(4)*zsize*3/4*pi) ' mm^3']);    
end
% classes of tissue
GM = 1; WM = 2; GMR = 3; GMSR = 4; WMR = 5; WMSR = 6; AI = 7; VO = 8; CSF = 9; HEGM = 10; HEWM = 11;
if handles.annotationMode ~= 3
    mask = h.createMask;
    delete(h);
    brainslice = squeeze(handles.brain(:,:,handles.sliceNum));
    brainslice_o = squeeze(handles.brain_o(:,:,handles.sliceNum));
    switch type
    case 0
        brainslice( (brainslice_o == GM) & mask) = GMR;
        brainslice( (brainslice_o == WM) & mask) = WMR;  
    case 1
        brainslice( (brainslice_o == GM) & mask) = GMSR;
        brainslice( (brainslice_o == WM) & mask) = WMSR;
    case 2
        brainslice( (brainslice == GM) & mask) = HEWM;
        brainslice( (brainslice == WM) & mask) = HEGM;  
    end
    handles.brain(:,:,handles.sliceNum) = brainslice;
else
    delete(h);
    mask = create3DEllipseMask(position(1), position(2), zpos-zsize/2, position(3), position(4), zsize);    
    switch type
    case 0
        handles.brain( (handles.brain_o == GM) & mask) = GMR;
        handles.brain( (handles.brain_o == WM) & mask) = WMR;       
    case 1    
        handles.brain( (handles.brain_o == GM) & mask) = GMSR;
        handles.brain( (handles.brain_o == WM) & mask) = WMSR;
    case 2
        handles.brain( (handles.brain == GM) & mask) = HEGM;
        handles.brain( (handles.brain == WM) & mask) = HEWM;        
    end    
end

function mask = create3DEllipseMask(xmin, ymin, zmin, xsize, ysize, zsize)
xmid = xmin+xsize/2;
ymid = ymin+ysize/2;
zmid = zmin+zsize/2;
xmin = floor(xmin);
ymin = floor(ymin);
zmin = floor(zmin);
if xmin < 1
    xmin = 1;
end
xmax = xmin+ceil(xsize);
if xmax > 256
    xmax = 256;
end
if ymin < 1
    ymin = 1;
end
ymax = ymin+ceil(ysize);
if ymax > 256
    ymax = 256;
end
if zmin < 1
    zmin = 1;
end
zmax = zmin+ceil(zsize);
if zmax > 256
    zmax = 256;
end
mask = zeros(256,256,256);
a = (xsize*xsize)/4;
b = (ysize*ysize)/4;
c = (zsize*zsize)/4;
for z=zmin:zmax
    for y=ymin:ymax
        for x=xmin:xmax
            v = (x-xmid)*(x-xmid)/a + (y-ymid)*(y-ymid)/b + (z-zmid)*(z-zmid)/c;
            mask(y,x,z) = (v  <= 1);
        end
    end
end

function [segbrain_o segbrain mrbrain] = loadBrain()
if ~exist('brain.raw','file')
    error('brain.raw not found!');
end
if ~exist('mrbrain.raw','file')
    error('mrbrain.raw not found!');
end

fid = fopen('brain.raw','rb');
segbrain_o = zeros(256,256,256);
for i=1:256
    segbrain_o(:,:,i) = fread(fid,[256 256],'char');
end
fclose(fid); 

if exist('brain.mat','file')
    s = load('brain.mat','brain');
    segbrain = s.brain;
else
    segbrain = segbrain_o;
end

fid = fopen('mrbrain.raw','rb');
mrbrain = zeros(256,256,256);
for i=1:256
    mrbrain(:,:,i) = fread(fid,[256 256],'uint16');
end
fclose(fid);
