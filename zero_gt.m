clear all;

% This script acts as a benchmark in which to compare the 4 approximations;
% it is a ground-truth realisation of the STF, modelling it as a 256-tap
% FIR at each point in space.

%|-----------Parameters----------|

c = 343; % speed of sound in metres per second
theta_s = -pi/2; % azimuth start-point of source
theta_e = pi/2; % azimuth end-point of source
Res = 100; % movement 'resolution' (number of times azimuth is varied, thus number of filter sets generated)
gamma_le = -pi/2; % azimuth of left ear
gamma_re = pi/2; % azimuth of right ear
e2e = 0.215; % distance between ears, used for ITDs
a_h = 0.0875; % average radius of adult human head (according to B & D)
rho_d = 32;
threshold = 0.001;
fir_order = 256;
f = linspace(1,22000,1000); % frequency discretization

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

%|----------Initilizations----------|

h = waitbar(0,'Initializing...');

r = rho_d*a_h;
theta = linspace(theta_s,theta_e,Res);
Nyq = Fs/2;
w0 = c/a_h;
gamma = theta;

y = zeros(L,2);
seg = floor(L/Res);

f_n = (f/max(f));
f_n(1) = 0;

tic;

for k=1:Res

%|----------AOIs and ITD calculations----------|  

if sign(gamma(k))==-1 % i.e. source is in front-left quarter space
    theta_l = -gamma_le-abs(gamma(k)); % the left ear is the ipsilateral ear
    theta_r = gamma_re+abs(gamma(k)); % the right ear is the contralateral ear
    Td_l = -(e2e/(2*c))*cos(theta_l)+(e2e/(2*c)); % delay between centre and left ear
    Td_r = (e2e/(2*c))*(abs(theta_r)-(pi/2))+(e2e/(2*c)); % delay between centre and right ear
else % i.e. source is in front-right quarter space
    theta_l = -gamma_le+abs(gamma(k));
    theta_r = gamma_re-abs(gamma(k)); % vice versa
    Td_l = (e2e/(2*c))*(abs(theta_l)-(pi/2))+(e2e/(2*c)); % time taken for signal to reach left ear
    Td_r = -(e2e/(2*c))*cos(theta_r)+(e2e/(2*c)); % time taken for signal to reach right ear
end

%|----------Spherical frequency response----------|

for i=1:length(f)
    H_STF_l(i) = sphere(a_h,r,theta_l,f(i),threshold);
    H_STF_r(i) = sphere(a_h,r,theta_r,f(i),threshold);
end

%|----------Head shadow IIR (Binaural)----------|

b_l = fir2(fir_order,f_n,abs(H_STF_l));
b_r = fir2(fir_order,f_n,abs(H_STF_r));

% HS_IIR_l = dfilt.df1(b_l,a_l);
% HS_IIR_r = dfilt.df1(b_r,a_r);

HS_FIR_l = dfilt.dffir(b_l);
HS_FIR_r = dfilt.dffir(b_r);

%|----------ITD----------|

Td_l_n = Td_l*Fs; % delay in 'floating sample points' (for left ear)
Td_r_n = Td_r*Fs; % delay in 'floating sample points' (for right ear)

Td_l_i = floor(Td_l_n); % integer portion of delay in samples (for left ear)
Td_r_i = floor(Td_r_n); % integer portion of delay in samples (for right ear)

Hd_l = dfilt.delay(Td_l_i);
Hd_r = dfilt.delay(Td_r_i);

H_l = dfilt.cascade(HS_FIR_l,Hd_l);
H_r = dfilt.cascade(HS_FIR_r,Hd_r);

%|----------Signal flow----------|

if k==1
    y(1:seg,1) = filter(H_l,x(1:seg)); % Headshadow and ITD
    y(1:seg,2) = filter(H_r,x(1:seg));
elseif k<Res
    y(k*seg:(k+1)*seg,1) = filter(H_l,x(k*seg:(k+1)*seg)); % Headshadow and ITD
    y(k*seg:(k+1)*seg,2) = filter(H_r,x(k*seg:(k+1)*seg));
end
if k==Res
    y(end-seg:end,1) = filter(H_l,x(end-seg:end)); % Headshadow and ITD
    y(end-seg:end,2) = filter(H_r,x(end-seg:end));
end

waitbar(k/Res,h,'Working...')

end

toc;

close(h)

% |----------WAV output----------|

rho_str = num2str(rho_d);

audiowrite(['zero_gt_',rho_str,'.wav'],y,Fs,'BitsPerSample',16);




