function varargout = GUIDE_HELP(varargin)
% GUIDE_HELP MATLAB code for GUIDE_HELP.fig
%      GUIDE_HELP, by itself, creates a new GUIDE_HELP or raises the existing
%      singleton*.
%
%      H = GUIDE_HELP returns the handle to a new GUIDE_HELP or the handle to
%      the existing singleton*.
%
%      GUIDE_HELP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIDE_HELP.M with the given input arguments.
%
%      GUIDE_HELP('Property','Value',...) creates a new GUIDE_HELP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUIDE_HELP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUIDE_HELP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUIDE_HELP

% Last Modified by GUIDE v2.5 05-Dec-2024 20:22:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUIDE_HELP_OpeningFcn, ...
                   'gui_OutputFcn',  @GUIDE_HELP_OutputFcn, ...
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


% --- Executes just before GUIDE_HELP is made visible.
function GUIDE_HELP_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to GUIDE_HELP (see VARARGIN)
    
    % Choose default command line output for GUIDE_HELP
    handles.output = hObject;
    
    % Update handles structure
    guidata(hObject, handles);
    
    % UIWAIT makes GUIDE_HELP wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

    thresholdingImg = imread('Symbols\thresholding.png');
    axes(handles.axes1);
    imshow(thresholdingImg);


% --- Outputs from this function are returned to the command line.
function varargout = GUIDE_HELP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
