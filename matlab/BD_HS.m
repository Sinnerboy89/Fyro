clear all;

%|-----------Parameters----------|

c = 343; % speed of sound in metres per second
theta = pi/4; % azimuth of source
gamma_le = -pi/2; % azimuth of left ear
gamma_re = pi/2; % azimuth of right ear
a_h = 0.0875; % average radius of adult human head (according to B & D)
e2e = 0.215; % distance between ears, used for ITDs
alpha_min = 0.1; % values chosen by B & D as reasonable starting points before personalization
theta_min = 150*pi/180; % values chosen by B & D as reasonable starting points before personalization

%|----------Audio file read-in----------|

currentFolder = pwd;

% Get the full path of the audio file that the user wants to binauralize

audiofilepath = fullfile(currentFolder, '*.*');
disp('Pick an audio sample please')
[audiofilename, folder] = uigetfile(audiofilepath, 'Pick a sample to binauralize');
if audiofilename == 0  % user clicked cancel
    error('Fine then. Dont use my script, see if I care.');
end
filepathin = fullfile(folder, audiofilename);

[x,Fs] = audioread(filepathin);

if size(x,2)>1
    x = (x(:,1)+x(:,2))/2;
    disp('input is stereo - converted to mono')
end

L = length(x);
y = zeros(L,2);

tic;

%|----------Initilizations----------|

Nyq = Fs/2;
w0 = c/a_h;
gamma = theta;

%|----------AOIs and ITD calculations----------|  

if sign(gamma)==-1 % i.e. source is in front-left quarter space
    theta_l = -gamma_le-abs(gamma); % the left ear is the ipsilateral ear
    theta_r = gamma_re+abs(gamma); % the right ear is the contralateral ear
    Td_l = -(e2e/(2*c))*cos(theta_l)+(e2e/(2*c)); % delay between centre and left ear
    Td_r = (e2e/(2*c))*(abs(theta_r)-(pi/2))+(e2e/(2*c)); % delay between centre and right ear
else % i.e. source is in front-right quarter space
    theta_l = -gamma_le+abs(gamma);
    theta_r = gamma_re-abs(gamma); % vice versa
    Td_l = (e2e/(2*c))*(abs(theta_l)-(pi/2))+(e2e/(2*c)); % time taken for signal to reach left ear
    Td_r = -(e2e/(2*c))*cos(theta_r)+(e2e/(2*c)); % time taken for signal to reach right ear
end

%|----------Spherical frequency response approximation----------|

%alpha represents an approximation to the frequency response of an ideal
%rigid sphere (as described by Rayleigh), and is valid for 0<abs(gamma)<5*pi/6. Beyond this angle, response
%starts to rise to the "bright spot" at pi - but this is non-linear, and not simulated here.

alpha_l = (1+alpha_min/2)+(1-alpha_min/2)*cos(theta_l*pi/theta_min);
alpha_r = (1+alpha_min/2)+(1-alpha_min/2)*cos(theta_r*pi/theta_min);

%|----------Head shadow IIR (Binaural)----------|

b_l = [w0+(alpha_l*Fs),w0-(alpha_l*Fs)];
b_r = [w0+(alpha_r*Fs),w0-(alpha_r*Fs)];
a = [w0+Fs,w0-Fs];

HS_IIR_l = dfilt.df1(b_l,a);
HS_IIR_r = dfilt.df1(b_r,a);

%|----------Signal flow----------|

y(:,1) = filter(HS_IIR_l,x); % Headshadow
y(:,2) = filter(HS_IIR_r,x);

toc;

soundsc(y,Fs)




