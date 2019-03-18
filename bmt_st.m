% clear all;

% BMT-based azimuth/range synthesis (stationary source)

% parameters
c = 343;
theta = -pi/2;
gamma_le = -pi/2; % azimuth location of left ear
gamma_re = pi/2; % azimuth location of right ear
e2e = 0.215; % distance between ears (used for ITDs)
rho_d = 1.25;

% IIR filterbank read-in
currentFolder = pwd;
iirfilepath = fullfile(currentFolder, '*.*');
disp('Pick IIR filterbank please')
[iirfilename, folder] = uigetfile(iirfilepath, 'Pick IIR filterbank please');
if iirfilename == 0
    error('Fine then. Dont use my script, see if I care.');
end
iirfilepathin = fullfile(folder, iirfilename);
load(iirfilepathin);

% audio file read-in (convert to mono if required)
audiofilepath = fullfile(currentFolder, '*.*');
disp('Pick an audio sample please')
[audiofilename, folder] = uigetfile(audiofilepath, 'Pick a sample to binauralize');
if audiofilename == 0
    error('Fine then. Dont use my script, see if I care.');
end
audiofilepathin = fullfile(folder, audiofilename);
[x,Fs] = audioread(audiofilepathin);
if size(x,2)>1
    x = (x(:,1)+x(:,2))/2;
    disp('input is stereo - converted to mono')
end
L = length(x);

% output signal init
y = zeros(L,2);

% AOI and ITD calculations
if sign(theta)==-1 % i.e. source is in front-left quarter space
    theta_l = -gamma_le-abs(theta); % the left ear is the ipsilateral ear
    theta_r = gamma_re+abs(theta); % the right ear is the contralateral ear
    Td_l = -(e2e/(2*c))*cos(theta_l)+(e2e/(2*c)); % delay between centre and left ear
    Td_r = (e2e/(2*c))*(abs(theta_r)-(pi/2))+(e2e/(2*c)); % delay between centre and right ear
else % i.e. source is in front-right quarter space
    theta_l = -gamma_le+abs(theta);
    theta_r = gamma_re-abs(theta); % vice versa
    Td_l = (e2e/(2*c))*(abs(theta_l)-(pi/2))+(e2e/(2*c)); % time taken for signal to reach left ear
    Td_r = -(e2e/(2*c))*cos(theta_r)+(e2e/(2*c)); % time taken for signal to reach right ear
end

% ITD as all-pass delay filter
Td_l_n = Td_l*Fs; % delay in 'floating sample points' (for left ear)
Td_r_n = Td_r*Fs; % delay in 'floating sample points' (for right ear)
Td_l_i = floor(Td_l_n); % integer portion of delay in samples (for left ear)
Td_r_i = floor(Td_r_n); % integer portion of delay in samples (for right ear)
Hd_l = dfilt.delay(Td_l_i);
Hd_r = dfilt.delay(Td_r_i);

% find 'nearest' filters for ILD (hardcode discretization for now)
azi_density = 91;
rho_density = 10;
theta_grid = linspace(0, pi, azi_density);
rho_grid = 1.15.^(ones(1, rho_density)+(linspace(0, 249, rho_density)/10));
[~, theta_l_idx] = min(abs(theta_grid - theta_l));
[~, theta_r_idx] = min(abs(theta_grid - theta_r));
[~, rho_idx] = min(abs(rho_grid - rho_d));
HS_IIR_l = H_IIR(theta_l_idx, rho_idx);
HS_IIR_r = H_IIR(theta_r_idx, rho_idx);

% cascading ILD IIR with ITD delay
HS_l = dfilt.cascade(HS_IIR_l,Hd_l);
HS_r = dfilt.cascade(HS_IIR_r,Hd_r);

% apply cascade to input signal
y(:,1) = filter(HS_l,x);
y(:,2) = filter(HS_r,x);

% normalize
norm = max(max(y));
y = y / norm;

% output
theta_str = num2str(theta);
rho_str = num2str(rho_d);
audiowrite(['bmt_a', theta_str, '_r', rho_str, '.wav'], y, Fs, 'BitsPerSample', 16);

