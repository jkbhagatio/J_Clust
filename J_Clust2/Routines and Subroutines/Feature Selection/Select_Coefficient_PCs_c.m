function varargout = Select_Coefficient_PCs_c(varargin)

%Select_Coefficient_PCs_c MATLAB code file for Select_Coefficient_PCs_c.fig
%
%Description: The call of this .m file, within the main J_Clust2 GUI, by selection of 'Concatenated PC Select' within the 
%'Feature Coefficient Selection'  button group, launches a pop-up sub-GUI, which allows the user to select which principal component to use,
%to calculate principal component scores of all concatenated spike waveforms corresponding to the selected principal component
%
% Begin initialization code - DO NOT EDIT
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Select_Coefficient_PCs_c_OpeningFcn, ...
                   'gui_OutputFcn',  @Select_Coefficient_PCs_c_OutputFcn, ...
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


% --- Executes just before Select_Coefficient_PCs_c is made visible.
function Select_Coefficient_PCs_c_OpeningFcn(select_coeffs_pc_c_obj, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% select_coeffs_pc_c_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Select_Coefficient_PCs_c (see VARARGIN)

% Choose default command line output for Select_Coefficient_PCs_c
handles.output = select_coeffs_pc_c_obj;

% Update handles structure
J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);

handles.PC_coeffs_c = data_from_main.features{11}; %PC Coeffs in format [num_samples X PC#]
handles.PC_var_c = data_from_main.features{12};

guidata(select_coeffs_pc_c_obj, handles);
set(select_coeffs_pc_c_obj,'menubar','figure')

% J_Clust_tag = findobj('Tag','J_Clust_fig');
% data_from_main = guidata(J_Clust_tag);
% 
% switch data_from_main.feat_disp1
%     case 'PCA Scores'
%         handles.coeffs = data_from_main.features{8};
%     case 'PCA Scores for Concatenated Waveforms'
%         handles.coeffs = data_from_main.features{10};
%     case 'Wavelet Coefficients'
%         %...
%     case 'Wavelet Coefficients for Concatenated Waveforms'
% end

% UIWAIT makes Select_Coefficient_PCs_c wait for user response (see UIRESUME)
% uiwait(handles.select_coeffs_pc_c_fig);


% --- Outputs from this function are returned to the command line.
function varargout = Select_Coefficient_PCs_c_OutputFcn(select_coeffs_pc_c_obj, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% select_coeffs_pc_c_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in selection_menu1.
function selection_menu1_Callback(select_coeffs_pc_c_obj, eventdata, handles)
% select_coeffs_pc_c_obj    handle to selection_menu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(select_coeffs_pc_c_obj,'String')) returns selection_menu1 contents as cell array
%        contents{get(select_coeffs_pc_c_obj,'Value')} returns selected item from selection_menu1

coeff_contents = get(select_coeffs_pc_c_obj, 'String');
cur_coeff = coeff_contents(get(select_coeffs_pc_c_obj, 'Value'));
if strcmp(cur_coeff, '') || strcmp(cur_coeff, 'Select Coefficient 1')
    return;
end
handles.cur_coeff1 = str2num(cell2mat(cur_coeff));

guidata(select_coeffs_pc_c_obj, handles);

axes(handles.coeff1_plot)
plot(handles.PC_coeffs_c(:,handles.cur_coeff1))
title(['Variance Accounted for: ', num2str(handles.PC_var_c(handles.cur_coeff1,1)), ' Percent'], 'FontSize', 8.5)



% --- Executes during object creation, after setting all properties.
function selection_menu1_CreateFcn(select_coeffs_pc_c_obj, eventdata, handles)
% select_coeffs_pc_c_obj    handle to selection_menu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(select_coeffs_pc_c_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(select_coeffs_pc_c_obj,'BackgroundColor','white');
end


% --- Executes on selection change in selection_menu2.
function selection_menu2_Callback(select_coeffs_pc_c_obj, eventdata, handles)
% select_coeffs_pc_c_obj    handle to selection_menu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(select_coeffs_pc_c_obj,'String')) returns selection_menu2 contents as cell array
%        contents{get(select_coeffs_pc_c_obj,'Value')} returns selected item from selection_menu2

coeff_contents = get(select_coeffs_pc_c_obj, 'String');
cur_coeff = coeff_contents(get(select_coeffs_pc_c_obj, 'Value'));
if strcmp(cur_coeff, '') || strcmp(cur_coeff, 'Select Coefficient 2')
    return;
end
handles.cur_coeff2 = str2num(cell2mat(cur_coeff));

guidata(select_coeffs_pc_c_obj, handles);

axes(handles.coeff2_plot)
plot(handles.PC_coeffs_c(:,handles.cur_coeff2))
title(['Variance Accounted for: ', num2str(handles.PC_var_c(handles.cur_coeff2,1)), ' Percent'], 'FontSize', 8.5)



% --- Executes during object creation, after setting all properties.
function selection_menu2_CreateFcn(select_coeffs_pc_c_obj, eventdata, handles)
% select_coeffs_pc_c_obj    handle to selection_menu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(select_coeffs_pc_c_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(select_coeffs_pc_c_obj,'BackgroundColor','white');
end


% --- Executes on selection change in selection_menu3.
function selection_menu3_Callback(select_coeffs_pc_c_obj, eventdata, handles)
% select_coeffs_pc_c_obj    handle to selection_menu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(select_coeffs_pc_c_obj,'String')) returns selection_menu3 contents as cell array
%        contents{get(select_coeffs_pc_c_obj,'Value')} returns selected item from selection_menu3

coeff_contents = get(select_coeffs_pc_c_obj, 'String');
cur_coeff = coeff_contents(get(select_coeffs_pc_c_obj, 'Value'));
if strcmp(cur_coeff, '') || strcmp(cur_coeff, 'Select Coefficient 3')
    return;
end
handles.cur_coeff3 = str2num(cell2mat(cur_coeff));

guidata(select_coeffs_pc_c_obj, handles);

axes(handles.coeff3_plot)
plot(handles.PC_coeffs_c(:,handles.cur_coeff3))
title(['Variance Accounted for: ', num2str(handles.PC_var_c(handles.cur_coeff3,1)), ' Percent'], 'FontSize', 8.5)


% --- Executes during object creation, after setting all properties.
function selection_menu3_CreateFcn(select_coeffs_pc_c_obj, eventdata, handles)
% select_coeffs_pc_c_obj    handle to selection_menu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(select_coeffs_pc_c_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(select_coeffs_pc_c_obj,'BackgroundColor','white');
end



% --- Executes on selection change in selection_menu4.
function selection_menu4_Callback(select_coeffs_pc_c_obj, eventdata, handles)
% select_coeffs_pc_c_obj    handle to selection_menu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(select_coeffs_pc_c_obj,'String')) returns selection_menu4 contents as cell array
%        contents{get(select_coeffs_pc_c_obj,'Value')} returns selected item from selection_menu4

coeff_contents = get(select_coeffs_pc_c_obj, 'String');
cur_coeff = coeff_contents(get(select_coeffs_pc_c_obj, 'Value'));
if strcmp(cur_coeff, '') || strcmp(cur_coeff, 'Select Coefficient 4')
    return;
end
handles.cur_coeff4 = str2num(cell2mat(cur_coeff)); 

guidata(select_coeffs_pc_c_obj, handles);

axes(handles.coeff4_plot)
plot(handles.PC_coeffs_c(:,handles.cur_coeff4))
title(['Variance Accounted for: ', num2str(handles.PC_var_c(handles.cur_coeff4,1)), ' Percent'], 'FontSize', 8.5)



% --- Executes during object creation, after setting all properties.
function selection_menu4_CreateFcn(select_coeffs_pc_c_obj, eventdata, handles)
% select_coeffs_pc_c_obj    handle to selection_menu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(select_coeffs_pc_c_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(select_coeffs_pc_c_obj,'BackgroundColor','white');
end



% --- Executes on button press in confirm_button.
function confirm_button_Callback(select_coeffs_pc_c_obj, eventdata, handles)
% select_coeffs_pc_c_obj    handle to confirm_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global cur_coeffs_pcs_c
cur_coeffs_pcs_c = [handles.cur_coeff1; handles.cur_coeff2; handles.cur_coeff3; handles.cur_coeff4];
guidata(select_coeffs_pc_c_obj, handles);
close(Select_Coefficient_PCs_c)
