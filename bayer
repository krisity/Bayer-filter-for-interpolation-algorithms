function varargout = MSVIVI(varargin)
% MSVIVI MATLAB code for MSVIVI.fig
%      MSVIVI, by itself, creates a new MSVIVI or raises the existing
%      singleton*.
%
%      H = MSVIVI returns the handle to a new MSVIVI or the handle to
%      the existing singleton*.
%
%      MSVIVI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MSVIVI.M with the given input arguments.
%
%      MSVIVI('Property','Value',...) creates a new MSVIVI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MSVIVI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MSVIVI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MSVIVI

% Last Modified by GUIDE v2.5 27-Nov-2018 16:39:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MSVIVI_OpeningFcn, ...
                   'gui_OutputFcn',  @MSVIVI_OutputFcn, ...
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


% --- Executes just before MSVIVI is made visible.
function MSVIVI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MSVIVI (see VARARGIN)

% Choose default command line output for MSVIVI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MSVIVI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MSVIVI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function imageName_Callback(hObject, eventdata, handles)
% hObject    handle to imageName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageName as text
%        str2double(get(hObject,'String')) returns contents of imageName as a double


% --- Executes during object creation, after setting all properties.
function imageName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------
%OPEN FILE BUTTON
%------------------------------
% --- Executes on button press in browseImageButton.
function browseImageButton_Callback(hObject, eventdata, handles)
% hObject    handle to browseImageButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.png';'*.gif';'*.jpg';'*.tif';'*.bmp'},'File Selector');
image = strcat(pathname, filename);
axes(handles.axes1);
imshow(image)
zoom on 
set(handles.imageName,'string',filename);


%------------------------------
%BILINEAR INTERPOLATION BUTTON
%------------------------------
% --- Executes on button press in bilinearInterpolationButton.
function bilinearInterpolationButton_Callback(hObject, eventdata, handles)
% hObject    handle to bilinearInterpolationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = getimage(handles.axes1);
J = bilinear (image);
axes(handles.axes2);
imshow(J)  

image2 = getimage(handles.axes2);

[p,m] = calculation(image,image2);
set(handles.PSNRtextField,'string',p);
set(handles.MSEtextField,'string',mean(m));
set(handles.mseValues,'string',m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% function out = bilinear(in) 
% returns the color interpolated image using Bilinear 
% interpolation algorithm. 
% 
% Assumptions : in has following color patterns 
% 
% ------------------> x 
% | G R G R ... 
% | B G B G ... 
% | G R G R ... 
% | B G B G ... 
% | . . . . . 
% | . . . . . 
% | . . . . . 
% | % V y 
%
%
% Input : 
% 
% in : original image matrix (mxnx3), m&n even 
% 
% Output : 
% 
% out : color interpolated image 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = bilinear(in)
m = size(in,1); n = size(in,2);
inR = in(:,:,1); inG = in(:,:,2); inB = in(:,:,3); 
out = in; outR = inR; outG = inG; outB = inB;

% R channel 
for i=1:2:m-1, 
    outR(i,3:2:n-1) = 1/2*(inR(i,2:2:n-2)+inR(i,4:2:n));
end

for i=2:2:m-2, 
    outR(i,2:2:n) = 1/2*(inR(i-1,2:2:n)+inR(i+1,2:2:n)); 
    outR(i,3:2:n-1) = 1/4*(inR(i-1,2:2:n-2)+inR(i-1,4:2:n)+inR(i+1,2:2:n-2)+inR(i+1,4:2:n)); 
end

outR = round(outR); 
ind = find(outR>255); 
outR(ind) = 255;

% B channel 
for i=2:2:m, 
    outB(i,2:2:n-2) = 1/2*(inB(i,1:2:n-3)+inB(i,3:2:n-1));
end

for i=3:2:m-1, 
    outB(i,1:2:n-1) = 1/2*(inB(i-1,1:2:n-1)+inB(i+1,1:2:n-1)); 
    outB(i,2:2:n-2) = 1/4*(inB(i-1,1:2:n-3)+inB(i-1,3:2:n-1)+inB(i+1,1:2:n-3)+inB(i+1,3:2:n-1)); 
end

outB = round(outB); 
ind = outB>255;
outB(ind) = 255;

% G channel 
for i=2:2:m-2, 
    outG(i,3:2:n-1) = 0.25*(inG(i,2:2:n-2)+inG(i,4:2:n)+inG(i-1,3:2:n-1)+inG(i+1,3:2:n-1)); 
end

for i=3:2:m-1,
    outG(i,2:2:n-2) = 0.25*(inG(i,1:2:n-3)+inG(i,3:2:n-1)+inG(i-1,2:2:n-2)+inG(i+1,2:2:n-1)); 
end

outG(1,2:2:n-2) = 1/3*(inG(1,1:2:n-3)+inG(3:2:n-1)+inG(2,2:2:n-2)); 
outG(1,n) = 1/2*(inG(1,n-1)+inG(2,n)); 
outG(3:2:m-1,n) = 1/3*(inG(2:2:m-2,n)+inG(4:2:m,n)+inG(3:2:m-1,n-1)); 
outG(2:2:m-2,1) = 1/3*(inG(1:2:m-3,1)+inG(3:2:m-1,1)+inG(2:2:m-2,2));
outG(m,1) = 1/2*(inG(m-1,1)+inG(m,2)); 
outG(m,3:2:n-1) = 1/3*(inG(m,2:2:n-2)+inG(m,4:2:n)+inG(m-1,3:2:n-1));

outG = round(outG);
ind = find(outG>255);
outG(ind) = 255;

out(:,:,1) = outR; 
out(:,:,2) = outG; 
out(:,:,3) = outB;

%------------------------------
%MEDIAN INTERPOLATION BUTTON
%------------------------------
% --- Executes on button press in medianInterpolationButton.
function medianInterpolationButton_Callback(hObject, eventdata, handles)
% hObject    handle to medianInterpolationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = getimage(handles.axes1);

% J = medianInterpolation (image);
% axes(handles.axes2);
% imshow(J)  
% 
% image2 = getimage(handles.axes2);
% 
% [p,m] = calculation(image,image2);
% set(handles.PSNRtextField,'string',p);
% set(handles.MSEtextField,'string',mean(m));
% set(handles.mseValues,'string',m);


% axes(handles.axes2);
% imshow(B)  

function out = medianInterpolation(in)
m = size(in,1); 
n = size(in,2);
inR = in(:,:,1); 
inG = in(:,:,2); 
inB = in(:,:,3); 
out = in; 
outR = inR; 
outG = inG; 
outB = inB;

% R channel 
outR = medfilt2(outR);

% B channel 

outB = medfilt2(outB);
% G channel 
outG = medfilt2(outG);

out(:,:,1) = outR; 
out(:,:,2) = outG; 
out(:,:,3) = outB;


%------------------------------
%CHANG INTERPOLATION BUTTON
%------------------------------
% --- Executes on button press in changInterpolationButton.
function changInterpolationButton_Callback(hObject, eventdata, handles)
% hObject    handle to changInterpolationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = getimage(handles.axes1);
J = chang(image);
axes(handles.axes2);
imshow(J)  

image2 = getimage(handles.axes2);

[p,m] = calculation(image,image2);
set(handles.PSNRtextField,'string',p);
set(handles.MSEtextField,'string',mean(m));
set(handles.mseValues,'string',m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function out = vargra(in) 
% returns the color interpolated image using 
% interpolation based on variable number of gradients 
% method 
% 
% Assumptions : in has following color patterns 
% 
% ------------------> x 
% | G R G R ... 
% | B G B G ... 
% | G R G R ... 
% | B G B G ... 
% | . . . . . 
% | . . . . . 
% | . . . . . 
% | 
% V y 
% 
% 
% Input : 
% 
% in : original image matrix (mxnx3), m&n even 
% 
% Output : 
% 
% out : color interpolated image %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = chang(in)
m = size(in,1); n = size(in,2);
inR = in(:,:,1); inG = in(:,:,2); inB = in(:,:,3); 
out = in; 
outR = inR; outG = inG; outB = inB;

k1= 1.5; k2 = 0.5;

% Estimate the missing color values at a non-green pixel 
% First : consider at red pixels 
for i=3:2:m-3; 
    for j=4:2:n-2; 
% form 8 gradients : N,E,S,W,NE,SE,NW,SW 
        gra_N = abs(inG(i-1,j)-inG(i+1,j))+abs(inR(i-2,j)-inR(i,j))+1/2*abs(inB(i-1,j-1)-inB(i+1,j-1))+1/2*abs(inB(i-1,j+1)-inB(i+1,j+1))+1/2*abs(inG(i-2,j-1)-inG(i,j-1))+1/2*abs(inG(i-2,j+1)-inG(i,j+1));
        gra_E = abs(inG(i,j+1)-inG(i,j-1))+abs(inR(i,j+2)-inR(i,j))+1/2*abs(inB(i-1,j+1)-inB(i-1,j-1))+1/2*abs(inB(i+1,j+1)-inB(i+1,j-1))+1/2*abs(inG(i-1,j+2)-inG(i-1,j))+1/2*abs(inG(i+1,j+2)-inG(i+1,j));
        gra_S = abs(inG(i+1,j)-inG(i-1,j))+abs(inR(i+2,j)-inR(i,j))+1/2*abs(inB(i+1,j+1)-inB(i-1,j+1))+1/2*abs(inB(i+1,j-1)-inB(i-1,j-1))+1/2*abs(inG(i+2,j+1)-inG(i,j+1))+1/2*abs(inG(i+2,j-1)-inG(i,j-1));
        gra_W = abs(inG(i,j-1)-inG(i,j+1))+abs(inR(i,j-2)-inR(i,j))+1/2*abs(inB(i+1,j-1)-inB(i+1,j+1))+1/2*abs(inB(i-1,j-1)-inB(i-1,j+1))+1/2*abs(inG(i+1,j-2)-inG(i+1,j))+1/2*abs(inG(i-1,j-2)-inG(i-1,j));
        gra_NE = abs(inB(i-1,j+1)-inB(i+1,j-1))+abs(inR(i-2,j+2)-inR(i,j))+1/2*abs(inG(i-1,j)-inG(i,j-1))+1/2*abs(inG(i,j+1)-inG(i+1,j))+1/2*abs(inG(i-2,j+1)-inG(i-1,j))+1/2*abs(inG(i-1,j+2)-inG(i,j+1));
        gra_SE = abs(inB(i+1,j+1)-inB(i-1,j-1))+abs(inR(i+2,j+2)-inR(i,j))+1/2*abs(inG(i,j+1)-inG(i-1,j))+1/2*abs(inG(i+1,j)-inG(i,j-1))+1/2*abs(inG(i+1,j+2)-inG(i,j+1))+1/2*abs(inG(i+2,j+1)-inG(i+1,j));
        gra_NW = abs(inB(i-1,j-1)-inB(i+1,j+1))+abs(inR(i-2,j-2)-inR(i,j))+1/2*abs(inG(i-1,j)-inG(i,j+1))+1/2*abs(inG(i,j-1)-inG(i+1,j))+1/2*abs(inG(i-2,j-1)-inG(i-1,j))+1/2*abs(inG(i-1,j-2)-inG(i,j-1));
        gra_SW = abs(inB(i+1,j-1)-inB(i-1,j+1))+abs(inR(i+2,j-2)-inR(i,j))+1/2*abs(inG(i+1,j)-inG(i,j+1))+1/2*abs(inG(i,j-1)-inG(i-1,j))+1/2*abs(inG(i+2,j-1)-inG(i+1,j))+1/2*abs(inG(i+1,j-2)-inG(i,j-1));

% determine thresholds
gra = [gra_N gra_E gra_S gra_W gra_NE gra_SE gra_NW gra_SW];
gramax = max(gra); 
gramin = min(gra); 
T = k1*gramin+k2*(gramax-gramin); 
ind = find(gra<T); 
Rave = zeros(1,8); 
Gave = zeros(1,8); 
Bave = zeros(1,8); 
for k=1:length(ind),
    switch ind(k) 
        case 1, 
            Rave(1) = 1/2*(inR(i,j)+inR(i-2,j));
            Gave(1) = inG(i-1,j); 
            Bave(1) = 1/2*(inB(i-1,j-1)+inB(i-1,j+1)); 
        case 2, 
            Rave(2) = 1/2*(inR(i,j)+inR(i,j+2)); 
            Gave(2) = inG(i,j+1); 
            Bave(2) = 1/2*(inB(i-1,j+1)+inB(i+1,j+1)); 
        case 3, 
            Rave(3) = 1/2*(inR(i,j)+inR(i+2,j)); 
            Gave(3) = inG(i+1,j); 
            Bave(3) = 1/2*(inB(i+1,j-1)+inB(i+1,j+1)); 
        case 4, 
            Rave(4) = 1/2*(inR(i,j)+inR(i,j-2)); 
            Gave(4) = inG(i,j-1); 
            Bave(4) = 1/2*(inB(i-1,j-1)+inB(i+1,j-1)); 
        case 5, 
            Rave(5) = 1/2*(inR(i,j)+inR(i-2,j+2)); 
            Gave(5) = 1/4*(inG(i,j+1)+inG(i-1,j+2)+inG(i-1,j)+inG(i-2,j+1)); 
            Bave(5) = inB(i-1,j+1); 
        case 6, 
            Rave(6) = 1/2*(inR(i,j)+inR(i+2,j+2)); 
            Gave(6) = 1/4*(inG(i,j+1)+inG(i+1,j+2)+inG(i+1,j)+inG(i+2,j+1)); 
            Bave(6) = inB(i+1,j+1); 
        case 7, 
            Rave(7) = 1/2*(inR(i,j)+inR(i-2,j-2));
            Gave(7) = 1/4*(inG(i,j-1)+inG(i-1,j-2)+inG(i-1,j)+inG(i-2,j-1)); 
            Bave(7) = inB(i-1,j-1); 
        case 8, 
            Rave(8) = 1/2*(inR(i,j)+inR(i+2,j-2)); 
            Gave(8) = 1/4*(inG(i,j-1)+inG(i+1,j-2)+inG(i+1,j)+inG(i+2,j-1)); 
            Bave(8) = inB(i+1,j-1); 
    end
end
Rsum = sum(Rave); 
Gsum = sum(Gave); 
Bsum = sum(Bave); 
outG(i,j) = inR(i,j)+(Gsum-Rsum)/length(ind); 
outB(i,j) = inR(i,j)+(Bsum-Rsum)/length(ind); 
    end
end

% Second : consider at blue pixels
for i=4:2:m-2; 
    for j=3:2:n-3; 
% form 8 gradients : N,E,S,W,NE,SE,NW,SW
        gra_N = abs(inG(i-1,j)-inG(i+1,j))+abs(inB(i-2,j)-inB(i,j))+1/2*abs(inR(i-1,j-1)-inR(i+1,j-1))+1/2*abs(inR(i-1,j+1)-inR(i+1,j+1))+1/2*abs(inG(i-2,j-1)-inG(i,j-1))+1/2*abs(inG(i-2,j+1)-inG(i,j+1));
        gra_E = abs(inG(i,j+1)-inG(i,j-1))+abs(inB(i,j+2)-inB(i,j))+1/2*abs(inR(i-1,j+1)-inR(i-1,j-1))+1/2*abs(inR(i+1,j+1)-inR(i+1,j-1))+1/2*abs(inG(i-1,j+2)-inG(i-1,j))+1/2*abs(inG(i+1,j+2)-inG(i+1,j));
        gra_S = abs(inG(i+1,j)-inG(i-1,j))+abs(inB(i+2,j)-inB(i,j))+1/2*abs(inR(i+1,j+1)-inR(i-1,j+1))+1/2*abs(inR(i+1,j-1)-inR(i-1,j-1))+1/2*abs(inG(i+2,j+1)-inG(i,j+1))+1/2*abs(inG(i+2,j-1)-inG(i,j-1));
        gra_W = abs(inG(i,j-1)-inG(i,j+1))+abs(inB(i,j-2)-inB(i,j))+1/2*abs(inR(i+1,j-1)-inR(i+1,j+1))+1/2*abs(inR(i-1,j-1)-inR(i-1,j+1))+1/2*abs(inG(i+1,j-2)-inG(i+1,j))+1/2*abs(inG(i-1,j-2)-inG(i-1,j));
        gra_NE = abs(inR(i-1,j+1)-inR(i+1,j-1))+abs(inB(i-2,j+2)-inB(i,j))+1/2*abs(inG(i-1,j)-inG(i,j-1))+1/2*abs(inG(i,j+1)-inG(i+1,j))+1/2*abs(inG(i-2,j+1)-inG(i-1,j))+1/2*abs(inG(i-1,j+2)-inG(i,j+1));
        gra_SE = abs(inR(i+1,j+1)-inR(i-1,j-1))+abs(inB(i+2,j+2)-inB(i,j))+1/2*abs(inG(i,j+1)-inG(i-1,j))+1/2*abs(inG(i+1,j)-inG(i,j-1))+1/2*abs(inG(i+1,j+2)-inG(i,j+1))+1/2*abs(inG(i+2,j+1)-inG(i+1,j));
        gra_NW = abs(inR(i-1,j-1)-inR(i+1,j+1))+abs(inB(i-2,j-2)-inB(i,j))+1/2*abs(inG(i-1,j)-inG(i,j+1))+1/2*abs(inG(i,j-1)-inG(i+1,j))+1/2*abs(inG(i-2,j-1)-inG(i-1,j))+1/2*abs(inG(i-1,j-2)-inG(i,j-1));
        gra_SW = abs(inR(i+1,j-1)-inR(i-1,j+1))+abs(inB(i+2,j-2)-inB(i,j))+1/2*abs(inG(i+1,j)-inG(i,j+1))+1/2*abs(inG(i,j-1)-inG(i-1,j))+1/2*abs(inG(i+2,j-1)-inG(i+1,j))+1/2*abs(inG(i+1,j-2)-inG(i,j-1));

% determine thresholds
gra = [gra_N gra_E gra_S gra_W gra_NE gra_SE gra_NW gra_SW]; 
gramax = max(gra); 
gramin = min(gra); 
T = k1*gramin+k2*(gramax-gramin); 
ind = find(gra<T); 
Rave = zeros(1,8); Gave = zeros(1,8); Bave = zeros(1,8);
for k=1:length(ind), 
    switch ind(k) 
        case 1, 
            Bave(1) = 1/2*(inB(i,j)+inB(i-2,j)); 
            Gave(1) = inG(i-1,j); 
            Rave(1) = 1/2*(inR(i-1,j-1)+inR(i-1,j+1)); 
        case 2, 
            Bave(2) = 1/2*(inB(i,j)+inB(i,j+2)); 
            Gave(2) = inG(i,j+1); 
            Rave(2) = 1/2*(inR(i-1,j+1)+inR(i+1,j+1)); 
        case 3, 
            Bave(3) = 1/2*(inB(i,j)+inB(i+2,j)); 
            Gave(3) = inG(i+1,j); 
            Rave(3) = 1/2*(inR(i+1,j-1)+inR(i+1,j+1)); 
        case 4, 
            Bave(4) = 1/2*(inB(i,j)+inB(i,j-2)); 
            Gave(4) = inG(i,j-1); 
            Rave(4) = 1/2*(inR(i-1,j-1)+inR(i+1,j-1)); 
        case 5, 
            Bave(5) = 1/2*(inB(i,j)+inB(i-2,j+2)); 
            Gave(5) = 1/4*(inG(i,j+1)+inG(i-1,j+2)+inG(i-1,j)+inG(i-2,j+1)); 
            Rave(5) = inR(i-1,j+1); 
        case 6, 
            Bave(6) = 1/2*(inB(i,j)+inB(i+2,j+2)); 
            Gave(6) = 1/4*(inG(i,j+1)+inG(i+1,j+2)+inG(i+1,j)+inG(i+2,j+1)); 
            Rave(6) = inR(i+1,j+1);
        case 7, 
            Bave(7) = 1/2*(inB(i,j)+inB(i-2,j-2)); 
            Gave(7) = 1/4*(inG(i,j-1)+inG(i-1,j-2)+inG(i-1,j)+inG(i-2,j-1)); 
            Rave(7) = inR(i-1,j-1); 
        case 8, 
            Bave(8) = 1/2*(inB(i,j)+inB(i+2,j-2)); 
            Gave(8) = 1/4*(inG(i,j-1)+inG(i+1,j-2)+inG(i+1,j)+inG(i+2,j-1));
            Rave(8) = inR(i+1,j-1); 
    end
end
Rsum = sum(Rave); 
Gsum = sum(Gave); 
Bsum = sum(Bave); 
outG(i,j) = inB(i,j)+(Gsum-Bsum)/length(ind); 
outR(i,j) = inB(i,j)+(Rsum-Bsum)/length(ind); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Estimating the missing color values at the green pixel location 
% First : consider those green pixels at upper-left 2x2 corner 
for i=3:2:m-3; 
    for j=3:2:n-3; 
% form 8 gradients : N,E,S,W,NE,SE,NW,SW
        gra_N = abs(inB(i-1,j)-inB(i+1,j))+abs(inG(i-2,j)-inG(i,j))+1/2*abs(inG(i-1,j-1)-inG(i+1,j-1))+1/2*abs(inG(i-1,j+1)-inG(i+1,j+1))+1/2*abs(inR(i-2,j-1)-inR(i,j-1))+1/2*abs(inR(i-2,j+1)-inR(i,j+1));
        gra_E = abs(inR(i,j+1)-inR(i,j-1))+abs(inG(i,j+2)-inG(i,j))+1/2*abs(inG(i-1,j+1)-inG(i-1,j-1))+1/2*abs(inG(i+1,j+1)-inG(i+1,j-1))+1/2*abs(inB(i-1,j+2)-inB(i-1,j))+1/2*abs(inB(i+1,j+2)-inB(i+1,j));
        gra_S = abs(inB(i+1,j)-inB(i-1,j))+abs(inG(i+2,j)-inG(i,j))+1/2*abs(inG(i+1,j+1)-inG(i-1,j+1))+1/2*abs(inG(i+1,j-1)-inG(i-1,j-1))+1/2*abs(inR(i+2,j+1)-inR(i,j+1))+1/2*abs(inR(i+2,j-1)-inR(i,j-1));
        gra_W = abs(inR(i,j-1)-inR(i,j+1))+abs(inG(i,j-2)-inG(i,j))+1/2*abs(inG(i+1,j-1)-inG(i+1,j+1))+1/2*abs(inG(i-1,j-1)-inG(i-1,j+1))+1/2*abs(inB(i+1,j-2)-inB(i+1,j))+1/2*abs(inB(i-1,j-2)-inB(i-1,j));
        gra_NE = abs(inG(i-1,j+1)-inG(i+1,j-1))+abs(inG(i-2,j+2)-inG(i,j))+abs(inR(i-2,j+1)-inR(i,j-1))+abs(inB(i-1,j+2)-inB(i+1,j));
        gra_SE = abs(inG(i+1,j+1)-inG(i-1,j-1))+abs(inG(i+2,j+2)-inG(i,j))+abs(inB(i+1,j+2)-inB(i-1,j))+abs(inR(i+2,j+1)-inR(i,j-1));
        gra_NW = abs(inG(i-1,j-1)-inG(i+1,j+1))+abs(inG(i-2,j-2)-inG(i,j))+abs(inR(i-2,j-1)-inR(i,j+1))+abs(inB(i-1,j-2)-inB(i+1,j));
        gra_SW = abs(inG(i+1,j-1)-inG(i-1,j+1))+abs(inG(i+2,j-2)-inG(i,j))+abs(inR(i+2,j-1)-inR(i,j+1))+abs(inB(i+1,j-2)-inG(i-1,j));

% determine thresholds
gra = [gra_N gra_E gra_S gra_W gra_NE gra_SE gra_NW gra_SW]; 
gramax = max(gra);
gramin = min(gra); 
T = k1*gramin+k2*(gramax-gramin); 
ind = find(gra<T); 
Rave = zeros(1,8); Gave = zeros(1,8); Bave = zeros(1,8);
        for k=1:length(ind), 
            switch ind(k) 
                case 1,
                    Gave(1) = 1/2*(inG(i,j)+inG(i-2,j)); 
                    Bave(1) = inB(i-1,j); 
                    Rave(1) = 1/4*(inR(i-2,j-1)+inR(i-2,j+1)+inR(i,j-1)+inR(i,j+1)); 
                case 2, 
                    Gave(2) = 1/2*(inG(i,j)+inG(i,j+2)); 
                    Rave(2) = inR(i,j+1); 
                    Bave(2) = 1/4*(inB(i-1,j)+inB(i+1,j)+inB(i-1,j+2)+inB(i+1,j+2)); 
                case 3, 
                    Gave(3) = 1/2*(inG(i,j)+inG(i+2,j)); 
                    Bave(3) = inB(i+1,j); 
                    Rave(3) = 1/4*(inR(i,j-1)+inR(i,j+1)+inR(i+2,j-1)+inR(i+2,j+1)); 
                case 4, 
                    Gave(4) = 1/2*(inG(i,j)+inG(i,j-2)); 
                    Rave(4) = inR(i,j-1); 
                    Bave(4) = 1/4*(inB(i-1,j-2)+inB(i-1,j)+inB(i+1,j-2)+inB(i+1,j)); 
                case 5, 
                    Rave(5) = 1/2*(inR(i-2,j+1)+inR(i,j+1)); 
                    Bave(5) = 1/2*(inB(i-1,j)+inB(i-1,j+2)); 
                    Gave(5) = inG(i-1,j+1);
                case 6,
                    Rave(6) = 1/2*(inR(i,j+1)+inR(i+2,j+1)); 
                    Bave(6) = 1/2*(inB(i+1,j)+inB(i+1,j+2)); 
                    Gave(6) = inG(i+1,j+1); 
                case 7, 
                    Rave(7) = 1/2*(inR(i,j-1)+inR(i-2,j-1)); 
                    Bave(7) = 1/2*(inB(i-1,j-2)+inB(i-1,j)); 
                    Gave(7) = inG(i-1,j-1); 
                case 8, 
                    Rave(8) = 1/2*(inR(i,j-1)+inR(i+2,j-1)); 
                    Bave(8) = 1/2*(inB(i+1,j-2)+inB(i+1,j)); 
                    Gave(8) = inG(i+1,j-1); 
            end
        end
        Rsum = sum(Rave); 
        Gsum = sum(Gave); 
        Bsum = sum(Bave); 
        outR(i,j) = inG(i,j)+(Rsum-Gsum)/length(ind); 
        outB(i,j) = inG(i,j)+(Bsum-Gsum)/length(ind); 
    end
end

% Second : consider those green pixels at the lower-right corner
for i=4:2:m-2; 
    for j=4:2:n-2; 
% form 8 gradients : N,E,S,W,NE,SE,NW,SW
        gra_N = abs(inR(i-1,j)-inR(i+1,j))+abs(inG(i-2,j)-inG(i,j))+1/2*abs(inG(i-1,j-1)-inG(i+1,j-1))+1/2*abs(inG(i-1,j+1)-inG(i+1,j+1))+1/2*abs(inB(i-2,j-1)-inB(i,j-1))+1/2*abs(inB(i-2,j+1)-inB(i,j+1));
        gra_E = abs(inB(i,j+1)-inB(i,j-1))+abs(inG(i,j+2)-inG(i,j))+1/2*abs(inG(i-1,j+1)-inG(i-1,j-1))+1/2*abs(inG(i+1,j+1)-inG(i+1,j-1))+1/2*abs(inR(i-1,j+2)-inR(i-1,j))+1/2*abs(inR(i+1,j+2)-inR(i+1,j));
        gra_S = abs(inR(i+1,j)-inR(i-1,j))+abs(inG(i+2,j)-inG(i,j))+1/2*abs(inG(i+1,j+1)-inG(i-1,j+1))+1/2*abs(inG(i+1,j-1)-inG(i-1,j-1))+1/2*abs(inB(i+2,j+1)-inB(i,j+1))+1/2*abs(inB(i+2,j-1)-inB(i,j-1));
        gra_W = abs(inB(i,j-1)-inB(i,j+1))+abs(inG(i,j-2)-inG(i,j))+1/2*abs(inG(i+1,j-1)-inG(i+1,j+1))+1/2*abs(inG(i-1,j-1)-inG(i-1,j+1))+1/2*abs(inR(i+1,j-2)-inR(i+1,j))+1/2*abs(inR(i-1,j-2)-inR(i-1,j));
        gra_NE = abs(inG(i-1,j+1)-inG(i+1,j-1))+abs(inG(i-2,j+2)-inG(i,j))+abs(inB(i-2,j+1)-inB(i,j-1))+abs(inR(i-1,j+2)-inR(i+1,j));
        gra_SE = abs(inG(i+1,j+1)-inG(i-1,j-1))+abs(inG(i+2,j+2)-inG(i,j))+abs(inR(i+1,j+2)-inR(i-1,j))+abs(inB(i+2,j+1)-inB(i,j-1));
        gra_NW = abs(inG(i-1,j-1)-inG(i+1,j+1))+abs(inG(i-2,j-2)-inG(i,j))+abs(inB(i-2,j-1)-inB(i,j+1))+abs(inR(i-1,j-2)-inR(i+1,j));
        gra_SW = abs(inG(i+1,j-1)-inG(i-1,j+1))+abs(inG(i+2,j-2)-inG(i,j))+abs(inB(i+2,j-1)-inB(i,j+1))+abs(inR(i+1,j-2)-inR(i-1,j));

% determine thresholds
        gra = [gra_N gra_E gra_S gra_W gra_NE gra_SE gra_NW gra_SW]; 
        gramax = max(gra); 
        gramin = min(gra); 
        T = k1*gramin+k2*(gramax-gramin); 
        ind = find(gra<T);
Rave = zeros(1,8); Gave = zeros(1,8); Bave = zeros(1,8); 
for k=1:length(ind), 
    switch ind(k) 
        case 1, 
            Gave(1) = 1/2*(inG(i,j)+inG(i-2,j)); 
            Rave(1) = inR(i-1,j); 
            Bave(1) = 1/4*(inB(i-2,j-1)+inB(i-2,j+1)+inB(i,j-1)+inB(i,j+1)); 
        case 2, 
            Gave(2) = 1/2*(inG(i,j)+inG(i,j+2));
            Bave(2) = inB(i,j+1);
            Rave(2) = 1/4*(inR(i-1,j)+inR(i+1,j)+inR(i-1,j+2)+inR(i+1,j+2));
        case 3, 
            Gave(3) = 1/2*(inG(i,j)+inG(i+2,j)); 
            Rave(3) = inR(i+1,j); 
            Bave(3) = 1/4*(inB(i,j-1)+inB(i,j+1)+inB(i+2,j-1)+inB(i+2,j+1)); 
        case 4, 
            Gave(4) = 1/2*(inG(i,j)+inG(i,j-2)); 
            Bave(4) = inB(i,j-1); 
            Rave(4) = 1/4*(inR(i-1,j-2)+inR(i-1,j)+inR(i+1,j-2)+inR(i+1,j)); 
        case 5, 
            Bave(5) = 1/2*(inB(i-2,j+1)+inB(i,j+1)); 
            Rave(5) = 1/2*(inR(i-1,j)+inR(i-1,j+2)); 
            Gave(5) = inG(i-1,j+1); 
        case 6, 
            Bave(6) = 1/2*(inB(i,j+1)+inB(i+2,j+1)); 
            Rave(6) = 1/2*(inR(i+1,j)+inR(i+1,j+2)); 
            Gave(6) = inG(i+1,j+1); 
        case 7, 
            Bave(7) = 1/2*(inB(i,j-1)+inB(i-2,j-1)); 
            Rave(7) = 1/2*(inR(i-1,j-2)+inR(i-1,j)); 
            Gave(7) = inG(i-1,j-1); 
        case 8, 
            Bave(8) = 1/2*(inB(i,j-1)+inB(i+2,j-1)); 
            Rave(8) = 1/2*(inR(i+1,j-2)+inR(i+1,j)); 
            Gave(8) = inG(i+1,j-1); 
    end
end
Rsum = sum(Rave); 
Gsum = sum(Gave); 
Bsum = sum(Bave); 
outR(i,j) = inG(i,j)+(Rsum-Gsum)/length(ind); 
outB(i,j) = inG(i,j)+(Bsum-Gsum)/length(ind);
    end
end

outG = round(outG);
ind = outG>255; 
outG(ind) = 255; 
ind = outG<0; 
outG(ind) = 0;

outR = round(outR); 
ind = outR>255; 
outR(ind) = 255; 
ind = outR<0; 
outR(ind) = 0;

outB = round(outB); 
ind = outB>255; 
outB(ind) = 255;
ind = outB<0; 
outB(ind) = 0;

out(:,:,1) = outR; 
out(:,:,2) = outG; 
out(:,:,3) = outB;

%------------------------------
%LAROCHE PRESCOT INTERPOLATION BUTTON
%------------------------------
% --- Executes on button press in larochePrescotInterpolationButton.
function larochePrescotInterpolationButton_Callback(hObject, eventdata, handles)
% hObject    handle to larochePrescotInterpolationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = getimage(handles.axes1);
J = es2(image);
axes(handles.axes2);
imshow(J)  

image2 = getimage(handles.axes2);

[p,m] = calculation(image,image2);
set(handles.PSNRtextField,'string',p);
set(handles.MSEtextField,'string',mean(m));
set(handles.mseValues,'string',m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% function out = es2(in) 
% returns the color interpolated image using edge sensing 
% interpolation I 
% 
% Assumptions : in has following color patterns 
% 
% ------------------> x 
% | G R G R ... 
% | B G B G ... 
% | G R G R ... 
% | B G B G ... 
% | . . . . .
% | . . . . .
% | . . . . . 
% | 
% V y
% 
% 
% Input : 
% 
% in : original image matrix (mxnx3), m&n even 
% % Output : 
%
% out : color interpolated image 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function out = es2(in)
m = size(in,1); n = size(in,2); 
inR = in(:,:,1); inG = in(:,:,2); inB = in(:,:,3);
out = in; 
outR = inR; outG = inG; outB = inB;

% G channel 
for i=4:2:m-2,
    for j=3:2:n-3, 
        delta_H = abs(1/2*(inB(i,j-2)+inB(i,j+2))-inB(i,j)); 
        delta_V = abs(1/2*(inB(i-2,j)+inB(i+2,j))-inB(i,j)); 
        if delta_H < delta_V, 
            outG(i,j) = 1/2*(inG(i,j-1)+inG(i,j+1)); 
        elseif delta_H > delta_V,
            outG(i,j) = 1/2*(inG(i-1,j)+inG(i+1,j)); 
        else outG(i,j) = 1/4*(inG(i,j-1)+inG(i,j+1)+inG(i-1,j)+inG(i+1,j)); 
        end
    end
end

for i=3:2:m-3, 
    for j=4:2:n-2, 
        delta_H = abs(1/2*(inR(i,j-2)+inR(i,j+2))-inR(i,j));
        delta_V = abs(1/2*(inR(i-2,j)+inR(i+2,j))-inR(i,j)); 
        if delta_H < delta_V, 
            outG(i,j) = 1/2*(inG(i,j-1)+inG(i,j+1)); 
        elseif delta_H > delta_V, 
            outG(i,j) = 1/2*(inG(i-1,j)+inG(i+1,j));
        else outG(i,j) = 1/4*(inG(i,j-1)+inG(i,j+1)+inG(i-1,j)+inG(i+1,j)); 
        end
    end
end

outG = round(outG); 
ind = outG>255; 
outG(ind) = 255;

% R channel 
for i=1:2:m-1, 
    outR(i,3:2:n-1) = outG(i,3:2:n-1) + 1/2*(inR(i,2:2:n-2)-outG(i,2:2:n-2)+inR(i,4:2:n)-outG(i,4:2:n)); 
end

for i=2:2:m-2, 
    outR(i,2:2:n) = outG(i,2:2:n)+1/2*(inR(i-1,2:2:n)-outG(i-1,2:2:n)+inR(i+1,2:2:n)-outG(i+1,2:2:n)); 
    outR(i,3:2:n-1) = outG(i,3:2:n-1)+1/4*(inR(i-1,2:2:n-2)-outG(i-1,2:2:n-1)+inR(i-1,4:2:n)-outG(i-1,4:2:n)+inR(i+1,2:2:n-2)-outG(i+1,2:2:n-2)+inR(i+1,4:2:n)-outG(i+1,4:2:n)); 
end

outR = round(outR); 
ind = outR>255; 
outR(ind) = 255; 
ind = outR<0; 
outR(ind) = 0;

% B channel 
for i=2:2:m,
    outB(i,2:2:n-2) = outG(i,2:2:n-2)+1/2*(inB(i,1:2:n-3)-outG(i,1:2:n-3)+inB(i,3:2:n-1)-outG(i,3:2:n-1)); 
end

for i=3:2:m-1, 
    outB(i,1:2:n-1) = outG(i,1:2:n-1)+1/2*(inB(i-1,1:2:n-1)-outG(i-1,1:2:n-1)+inB(i+1,1:2:n-1)-outG(i+1,1:2:n-1)); 
    outB(i,2:2:n-2) = outG(i,2:2:n-2)+1/4*(inB(i-1,1:2:n-3)-outG(i-1,1:2:n-3)+inB(i-1,3:2:n-1)-outG(i-1,3:2:n-1)+inB(i+1,1:2:n-3)-outG(i+1,1:2:n-3)+inB(i+1,3:2:n-1)-outG(i+1,3:2:n-1)); 
end

outB = round(outB);
ind = outB>255; 
outB(ind) = 255;
ind = outB<0; 
outB(ind) = 0;

out(:,:,1) = outR; 
out(:,:,2) = outG; 
out(:,:,3) = outB;

%------------------------------
%HAMILTON ADAMS INTERPOLATION BUTTON
%------------------------------
% --- Executes on button press in HamiltonAdamsInterpolationButton.
function HamiltonAdamsInterpolationButton_Callback(hObject, eventdata, handles)
% hObject    handle to HamiltonAdamsInterpolationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = getimage(handles.axes1);
J = lcc1(image);
axes(handles.axes2);
imshow(uint8(J),'InitialMagnification',100);
image2 = getimage(handles.axes2);

[p,m] = calculation(image,image2);
set(handles.PSNRtextField,'string',p);
set(handles.MSEtextField,'string',mean(m));
set(handles.mseValues,'string',m);

% mse1 = mse(image,image2);
% set(handles.MSEtextField,'string',mse1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function out = lcc1(in) 
% returns the color interpolated image using Laplacian 
% second-order color correction I 
% 
% Assumptions : in has following color patterns 
% 
% ------------------> x 
% | G R G R ... 
% | B G B G ... 
% | G R G R ... 
% | B G B G ... 
% | . . . . . 
% | . . . . . 
% | . . . . . 
% | 
% V y 
%
% 
% Input : 
% 
% in : original image matrix (mxnx3), m&n even 
% 
% Output : 
% 
% out : color interpolated image 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = lcc1(in)
m = size(in,1); n = size(in,2); 
inR = in(:,:,1); inG = in(:,:,2); inB = in(:,:,3); 
out = in; 
outR = inR; outG = inG; outB = inB;

% G channel 
for i=4:2:m-2, 
    for j=3:2:n-3, 
        delta_H = abs(inB(i,j-2)+inB(i,j+2)-2*inB(i,j))+abs(inG(i,j-1)-inG(i,j+1)); 
        delta_V = abs(inB(i-2,j)+inB(i+2,j)-2*inB(i,j))+abs(inG(i-1,j)-inG(i+1,j)); 
        if delta_H < delta_V, 
            outG(i,j) = 1/2*(inG(i,j-1)+inG(i,j+1))+1/4*(2*inB(i,j)-inB(i,j-2)-inB(i,j+2)); 
        elseif delta_H > delta_V, 
            outG(i,j) = 1/2*(inG(i-1,j)+inG(i+1,j))+1/4*(2*inB(i,j)-inB(i-2,j)-inB(i+2,j)); 
        else outG(i,j) = 1/4*(inG(i,j-1)+inG(i,j+1)+inG(i-1,j)+inG(i+1,j))+1/8*(4*inB(i,j)-inB(i,j-2)-inB(i,j+2)-inB(i-2,j)-inB(i+2,j)); 
        end
    end
end

for i=3:2:m-3, 
    for j=4:2:n-2, 
        delta_H = abs(inR(i,j-2)+inR(i,j+2)-2*inR(i,j))+abs(inG(i,j-1)-inG(i,j+1)); 
        delta_V = abs(inR(i-2,j)+inR(i+2,j)-2*inR(i,j))+abs(inG(i-1,j)-inG(i+1,j));
        if delta_H < delta_V, 
            outG(i,j) = 1/2*(inG(i,j-1)+inG(i,j+1))+1/4*(2*inR(i,j)-inR(i,j-2)-inR(i,j+2)); 
        elseif delta_H > delta_V, 
            outG(i,j) = 1/2*(inG(i-1,j)+inG(i+1,j))+1/4*(2*inR(i,j)-inR(i-2,j)-inR(i+2,j)); 
        else
            outG(i,j) = 1/4*(inG(i,j-1)+inG(i,j+1)+inG(i-1,j)+inG(i+1,j))+1/8*(4*inR(i,j)-inR(i,j-2)-inR(i,j+2)-inR(i-2,j)-inR(i+2,j));
        end
    end
end

outG = round(outG); 
ind = find(outG>255); 
outG(ind) = 255;
ind = find(outG<0); 
outG(ind) = 0;

% R channel 
for i=1:2:m-1, 
    outR(i,3:2:n-1) = 1/2*(inR(i,2:2:n-2)+inR(i,4:2:n))+1/4*(2*outG(i,3:2:n-1)-outG(i,2:2:n-2)-outG(i,4:2:n)); 
end

for i=2:2:m-2, 
    outR(i,2:2:n) = 1/2*(inR(i-1,2:2:n)+inR(i+1,2:2:n))+1/4*(2*outG(i,2:2:n)-outG(i-1,2:2:n)-outG(i+1,2:2:n)); 
end

for i=2:2:m-2, 
    for j=3:2:n-1, 
        delta_P = abs(inR(i-1,j+1)-inR(i+1,j-1))+abs(2*outG(i,j)-outG(i-1,j+1)-outG(i+1,j-1)); 
        delta_N = abs(inR(i-1,j-1)-inR(i+1,j+1))+abs(2*outG(i,j)-outG(i-1,j-1)-outG(i+1,j+1)); 
        if delta_N < delta_P, 
            outR(i,j) = 1/2*(inR(i-1,j-1)+inR(i+1,j+1))+1/2*(2*outG(i,j)-outG(i-1,j-1)-outG(i+1,j+1)); 
        elseif delta_N > delta_P,
            outR(i,j) = 1/2*(inR(i-1,j+1)+inR(i+1,j-1))+1/2*(2*outG(i,j)-outG(i-1,j+1)-outG(i+1,j-1)); 
        else outR(i,j) = 1/4*(inR(i-1,j-1)+inR(i-1,j+1)+inR(i+1,j-1)+inR(i+1,j+1))+1/4*(4*outG(i,j)-outG(i-1,j-1)-outG(i-1,j+1)-outG(i+1,j-1)-outG(i+1,j+1)); 
        end
    end
end

outR = round(outR); 
ind = find(outR>255); 
outR(ind) = 255; 
ind = find(outR<0);
outR(ind) = 0;

% B channel 
for i=2:2:m, 
    outB(i,2:2:n-2) = 1/2*(inB(i,1:2:n-3)+inB(i,3:2:n-1))+1/4*(2*outG(i,2:2:n-2)-outG(i,1:2:n-3)-outG(i,3:2:n-1)); 
end

for i=3:2:m-1, 
    outB(i,1:2:n-1) = 1/2*(inB(i-1,1:2:n-1)+inB(i+1,1:2:n-1))+1/4*(2*outG(i,1:2:n-1)-outG(i-1,1:2:n-1)-outG(i+1,1:2:n-1)); 
end

for i=3:2:m-1,
    for j=2:2:n-2, 
        delta_P = abs(inB(i-1,j+1)-inB(i+1,j-1))+abs(2*outG(i,j)-outG(i-1,j+1)-outG(i+1,j-1)); 
        delta_N = abs(inB(i-1,j-1)-inB(i+1,j+1))+abs(2*outG(i,j)-outG(i-1,j-1)-outG(i+1,j+1)); 
        if delta_N < delta_P, 
            outB(i,j) = 1/2*(inB(i-1,j-1)+inB(i+1,j+1))+1/2*(2*outG(i,j)-outG(i-1,j-1)-outG(i+1,j+1)); 
        elseif delta_N > delta_P, 
            outB(i,j) = 1/2*(inB(i-1,j+1)+inB(i+1,j-1))+1/2*(2*outG(i,j)-outG(i-1,j+1)-outG(i+1,j-1)); 
        else outB(i,j) = 1/4*(inB(i-1,j-1)+inB(i-1,j+1)+inB(i+1,j-1)+inB(i+1,j+1))+1/4*(4*outG(i,j)-outG(i-1,j-1)-outG(i-1,j+1)-outG(i+1,j-1)-outG(i+1,j+1)); 
        end
    end
end

outB = round(outB);
ind = find(outB>255); 
outB(ind) = 255;
ind = find(outB<0);
outB(ind) = 0;

out(:,:,1) = outR; 
out(:,:,2) = outG; 
out(:,:,3) = outB;

function MSEtextField_Callback(hObject, eventdata, handles)
% hObject    handle to MSEtextField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MSEtextField as text
%        str2double(get(hObject,'String')) returns contents of MSEtextField as a double


% --- Executes during object creation, after setting all properties.
function MSEtextField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MSEtextField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PSNRtextField_Callback(hObject, eventdata, handles)
% hObject    handle to PSNRtextField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PSNRtextField as text
%        str2double(get(hObject,'String')) returns contents of PSNRtextField as a double


% --- Executes during object creation, after setting all properties.
function PSNRtextField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PSNRtextField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------
%SAVE BUTTON
%------------------------------
% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = getimage(handles.axes2);
filter = {'*.bmp';'*.tif'};
[file, path] = uiputfile(filter,'Save image as','interpolatedImage');
file=fullfile(path,file);
imwrite(image,file);

%------------------------------
%CLEAR BUTTON
%------------------------------
% --- Executes on button press in clearButton.
function clearButton_Callback(hObject, eventdata, handles)
% hObject    handle to clearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1,'reset')
cla(handles.axes2,'reset')
set(handles.PSNRtextField, 'String', '');
set(handles.MSEtextField, 'String', '');
set(handles.mseValues, 'string', ' ');
set(handles.imageName, 'String', '');

%------------------------------
%MSE AND PSNR CALCULATION
%------------------------------
function [peak, meanValue] = calculation(origImage, interpImage)

%Find MSE
[m,n] = size(origImage);
mse = (1/(m*n))*sum(sum((origImage-interpImage).^2));

%PSNR
mse1 = mean(mean((im2double(origImage) - im2double(interpImage)).^2, 1), 2);
psnr = 10 * log10(1 ./ mean(mse1,3));

peak = psnr;
meanValue = mse; 

fprintf('\nPSNR: %7.2f ', psnr);

%------------------------------
%PLOT BUTTON
%------------------------------
% --- Executes on button press in plotButton.
function plotButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
originalImage = getimage(handles.axes1);
[m,n] = size(originalImage);
bil = bilinear (originalImage);
ch = chang(originalImage);
lp = es2(originalImage);
ha = lcc1(originalImage);

[psnrBil, mseBil] = calculation(originalImage, bil);
[psnrCh, mseCh] = calculation(originalImage, ch);
[psnrLp, mseLp] = calculation(originalImage, lp);
[psnrHA, mseHA] = calculation(originalImage, ha);

% y = [psnrBil, psnrCh, psnrLp, psnrHA];
% y = [15.4788, 23.3839, 19.0141, 15.787];

% bilMse = [mseBil(:,:,1), mseBil(:,:,2), mseBil(:,:,3)];
% chMse = [mseCh(:,:,1), mseCh(:,:,2), mseCh(:,:,3)];
% lpMse = [mseLp(:,:,1), mseLp(:,:,2), mseLp(:,:,3)];
% haMse = [mseHA(:,:,1), mseHA(:,:,2), mseHA(:,:,3)];

bilMse = [mseBil(:,:,3), mseBil(:,:,2), mseBil(:,:,1)];
chMse = [mseCh(:,:,3), mseCh(:,:,2), mseCh(:,:,1)];
lpMse = [mseLp(:,:,3), mseLp(:,:,2), mseLp(:,:,1)];
haMse = [mseHA(:,:,3), mseHA(:,:,2), mseHA(:,:,1)];

figure
y = [bilMse; chMse; lpMse; haMse];
namesMse = {'Bilinear';'Chang';'Laroche-Prescot';'Hamilton-Adams'};
bar(y, 'group'); 
title('Bar Graph for individual MSE values');
legend('show');
legend('B','G','R')
set(gca,'xticklabel',namesMse)

figure
subplot(1,2,1);
namesPSNR = {'Bilinear';'Chang';'Laroche-Prescot';'Hamilton-Adams'};
x = [1:4]; 
z = [psnrBil, psnrCh, psnrLp, psnrHA];
bar(x,z);
text(1:length(z),z,num2str(z'),'vert','bottom','horiz','center'); 
title('Bar Graph for individual PSNR values');
set(gca,'xticklabel',namesPSNR);

subplot(1,2,2);
p = [mean(mseBil), mean(mseCh), mean(mseLp), mean(mseHA)];
nameMSE= {'Bilinear';'Chang';'Laroche-Prescot';'Hamilton-Adams'};
x = [1:4]; 
bar(x,p);
text(1:length(p),p,num2str(p'),'vert','bottom','horiz','center'); 
title('Bar Graph for MSE');
set(gca,'xticklabel',nameMSE);
