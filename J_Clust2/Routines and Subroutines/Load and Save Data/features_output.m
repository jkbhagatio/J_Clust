function varargout = features_output(varargin)

% features_output MATLAB code file for features_output.fig
%
%Description: The call of this .m file, within the main J_Clust2 GUI, by selection of 'All Unit Spike Features' within the 'Save Options' button group,
%launches a pop-up check-box menu, which allows the user to select which spike features to save for the units sorted in the current session
%
% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @features_output_OpeningFcn, ...
                   'gui_OutputFcn',  @features_output_OutputFcn, ...
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


% --- Executes just before features_output is made visible.
function features_output_OpeningFcn(features_output_obj, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% features_output_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to features_output (see VARARGIN)

% Choose default command line output for features_output
handles.output = features_output_obj;

% Update handles structure
handles.peak_amps_out = [];
handles.peak_peak_amps_out = [];
handles.power_out = [];
handles.pc_scores_out = [];
handles.pc_c_scores_out = [];
handles.wavelets_out = [];
handles.wavelets_c_out = [];
handles.width_out = [];


guidata(features_output_obj, handles);

% UIWAIT makes features_output wait for user response (see UIRESUME)
% uiwait(handles.features_output);


% --- Outputs from this function are returned to the command line.
function varargout = features_output_OutputFcn(features_output_obj, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% features_output_obj    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in peak_amps_output.
function peak_amps_output_Callback(features_output_obj, eventdata, handles)
% features_output_obj    handle to peak_amps_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(features_output_obj,'Value') returns toggle state of peak_amps_output

handles.peak_amps_out = get(features_output_obj,'Value');

guidata(features_output_obj, handles);


% --- Executes on button press in peak_peak_amps_output.
function peak_peak_amps_output_Callback(features_output_obj, eventdata, handles)
% features_output_obj    handle to peak_peak_amps_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(features_output_obj,'Value') returns toggle state of peak_peak_amps_output

handles.peak_peak_amps_out = get(features_output_obj,'Value');

guidata(features_output_obj, handles);


% --- Executes on button press in power_output.
function power_output_Callback(features_output_obj, eventdata, handles)
% features_output_obj    handle to power_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(features_output_obj,'Value') returns toggle state of power_output

handles.power_out = get(features_output_obj,'Value');

guidata(features_output_obj, handles);


% --- Executes on button press in pc_scores_output.
function pc_scores_output_Callback(features_output_obj, eventdata, handles)
% features_output_obj    handle to pc_scores_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(features_output_obj,'Value') returns toggle state of pc_scores_output

handles.pc_scores_out = get(features_output_obj,'Value');

guidata(features_output_obj, handles);


% --- Executes on button press in wavelet_coeffs_output.
function wavelet_coeffs_output_Callback(features_output_obj, eventdata, handles)
% features_output_obj    handle to wavelet_coeffs_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(features_output_obj,'Value') returns toggle state of wavelet_coeffs_output

handles.wavelets_out = get(features_output_obj,'Value');

guidata(features_output_obj, handles);


% --- Executes on button press in pc_c_scores_output.
function pc_c_scores_output_Callback(features_output_obj, eventdata, handles)
% features_output_obj    handle to pc_c_scores_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(features_output_obj,'Value') returns toggle state of pc_c_scores_output

handles.pc_c_scores_out = get(features_output_obj,'Value');

guidata(features_output_obj, handles);


% --- Executes on button press in wavelet_c_coeffs_output.
function wavelet_c_coeffs_output_Callback(features_output_obj, eventdata, handles)
% features_output_obj    handle to wavelet_c_coeffs_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(features_output_obj,'Value') returns toggle state of wavelet_c_coeffs_output

handles.wavelets_c_out = get(features_output_obj,'Value');

guidata(features_output_obj, handles);

function width_output_Callback(features_output_obj, eventdata, handles)
% features_output_obj    handle to width_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in width_output.
% Hint: get(features_output_obj,'Value') returns toggle state of width_output

handles.width_out = get(features_output_obj, 'Value');

guidata(features_output_obj, handles);




% --- Executes on button press in confirm.
function confirm_Callback(features_output_obj, eventdata, handles)
% features_output_obj    handle to confirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global features_vector_out
features_vector_out = zeros(8,1);

if handles.peak_amps_out
    features_vector_out(1) = 1;
end

if handles.peak_peak_amps_out
    features_vector_out(2) = 1;
end

if handles.power_out
    features_vector_out(3) = 1;
end

if handles.pc_scores_out
    features_vector_out(4) = 1;
end

if handles.pc_c_scores_out
    features_vector_out(5) = 1;
end

if handles.wavelets_out
    features_vector_out(6) = 1;
end

if handles.wavelets_c_out
    features_vector_out(7) = 1;
end

if handles.width_out
    features_vector_out(8) = 1;
end


close(features_output);
