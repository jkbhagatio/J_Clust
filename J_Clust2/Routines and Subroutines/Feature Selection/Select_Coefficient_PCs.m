function varargout = Select_Coefficient_PCs(varargin)

%Select_Coefficient_PCs MATLAB code file for Select_Coefficient_PCs.fig
%
%Description: The call of this .m file, within the main J_Clust2 GUI, by selection of 'PC Select' within the 'Feature Coefficient Selection' 
%button group, launches a pop-up sub-GUI, which allows the user to select which principal component to use, to calculate principal component scores
%of all spikes for all channels, corresponding to the selected principal component
%
% Begin initialization code - DO NOT EDIT
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Select_Coefficient_PCs_OpeningFcn, ...
                   'gui_OutputFcn',  @Select_Coefficient_PCs_OutputFcn, ...
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


% --- Executes just before Select_Coefficient_PCs is made visible.
function Select_Coefficient_PCs_OpeningFcn(select_coeffs_pc_obj, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% select_coeffs_pc_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Select_Coefficient_PCs (see VARARGIN)

% Choose default command line output for Select_Coefficient_PCs
handles.output = select_coeffs_pc_obj;

% Update handles structure
guidata(select_coeffs_pc_obj, handles);
set(select_coeffs_pc_obj,'menubar','figure')

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

% UIWAIT makes Select_Coefficient_PCs wait for user response (see UIRESUME)
% uiwait(handles.select_coeffs_pc_fig);


% --- Outputs from this function are returned to the command line.
function varargout = Select_Coefficient_PCs_OutputFcn(select_coeffs_pc_obj, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% select_coeffs_pc_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in selection_menu.
function selection_menu_Callback(select_coeffs_pc_obj, eventdata, handles)
% select_coeffs_pc_obj    handle to selection_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(select_coeffs_pc_obj,'String')) returns selection_menu contents as cell array
%        contents{get(select_coeffs_pc_obj,'Value')} returns selected item from selection_menu

%pass data from J_Clust to Select_Coefficient_PCs

J_Clust_tag = findobj('Tag','J_Clust_fig');
data_from_main = guidata(J_Clust_tag);

PC_coeffs = data_from_main.features{9}; %PC Coeffs in format [num_channels X num_samples X PC#]
PC_var = data_from_main.features{10};

coeff_contents = get(select_coeffs_pc_obj, 'String');
cur_coeff = coeff_contents(get(select_coeffs_pc_obj, 'Value'));
if strcmp(cur_coeff, '') || strcmp(cur_coeff, 'Select Coefficient')
    return;
end
handles.cur_coeff = str2num(cell2mat(cur_coeff));

guidata(select_coeffs_pc_obj, handles);

axes(handles.ch1_plot)
plot(PC_coeffs(1,:,handles.cur_coeff))
title(['Variance Accounted for: ', num2str(PC_var(handles.cur_coeff,1)), ' Percent'], 'FontSize', 8.5)
xlabel('Channel 1')

axes(handles.ch2_plot)
plot(PC_coeffs(2,:,handles.cur_coeff))
title(['Variance Accounted for: ', num2str(PC_var(handles.cur_coeff,2)), ' Percent'], 'FontSize', 8.5)
xlabel('Channel 2')

axes(handles.ch3_plot)
plot(PC_coeffs(3,:,handles.cur_coeff))
title(['Variance Accounted for: ', num2str(PC_var(handles.cur_coeff,3)), ' Percent'], 'FontSize', 8.5)
xlabel('Channel 3')

axes(handles.ch4_plot)
plot(PC_coeffs(4,:,handles.cur_coeff))
title(['Variance Accounted for: ', num2str(PC_var(handles.cur_coeff,4)), ' Percent'], 'FontSize', 8.5)
xlabel('Channel 4')
        
        %figure, subplot(1,2,2) first PC
        %pull-down menu to select up to 10th PC
        %title with variance contained for each PC and up to 10th
        %push button to select which PC to visualize in J_Clust
        %...


% --- Executes during object creation, after setting all properties.
function selection_menu_CreateFcn(select_coeffs_pc_obj, eventdata, handles)
% select_coeffs_pc_obj    handle to selection_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(select_coeffs_pc_obj,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(select_coeffs_pc_obj,'BackgroundColor','white');
end

% --- Executes on button press in confirm_button.
function confirm_button_Callback(select_coeffs_pc_obj, eventdata, handles)
% select_coeffs_pc_obj    handle to confirm_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global cur_coeff_pcs
cur_coeff_pcs = handles.cur_coeff;
guidata(select_coeffs_pc_obj, handles);
close(Select_Coefficient_PCs)
        
