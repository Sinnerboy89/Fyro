% clear all;
% 
% % BMT-based azimuth/range synthesis (linear azimuth source)
% 
% % parameters
% c = 343;
% theta_range = [-pi/2 pi/2];
% gamma_le = -pi/2; % azimuth location of left ear
% gamma_re = pi/2; % azimuth location of right ear
% e2e = 0.215; % distance between ears (used for ITDs)
% rho_d = 1.25;
puf = 0.05; % positional update frequency
% 
% % IIR filterbank read-in
% currentFolder = pwd;
% iirfilepath = fullfile(currentFolder, '*.*');
% disp('Pick IIR filterbank please')
% [iirfilename, folder] = uigetfile(iirfilepath, 'Pick IIR filterbank please');
% if iirfilename == 0
%     error('Fine then. Dont use my script, see if I care.');
% end
% iirfilepathin = fullfile(folder, iirfilename);
% load(iirfilepathin);
% 
% % audio file read-in (convert to mono if required)
% audiofilepath = fullfile(currentFolder, '*.*');
% disp('Pick an audio sample please')
% [audiofilename, folder] = uigetfile(audiofilepath, 'Pick a sample to binauralize');
% if audiofilename == 0
%     error('Fine then. Dont use my script, see if I care.');
% end
% audiofilepathin = fullfile(folder, audiofilename);
% [x_in,Fs] = audioread(audiofilepathin);
% if size(x_in,2)>1
%     x_in = (x_in(:,1)+x_in(:,2))/2;
%     disp('input is stereo - converted to mono')
% end
% L = length(x_in);

% output signal init
y_out = zeros(L,2);

% COLA setup
M = (Fs*puf) + 1;
R = (M-1)/2;
w = hann(M, 'periodic');
theta_ls = linspace(theta_range(1), theta_range(2), L);

for so=0:R:L-M
    % current signal window location
    ndx = so+1:so+M;
    x = x_in(ndx);
    
    % current azimuth
    theta = theta_ls(round(so + (M/2)));
    
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
    y_l = filter(HS_l,x);
    y_r = filter(HS_r,x);
    
    % window overlap-add
    y_out(ndx, 1) = y_out(ndx, 1) + y_l.*w;
    y_out(ndx, 2) = y_out(ndx, 2) + y_r.*w;
end

% normalize
norm = max(max(y_out));
y_out = y_out / norm;

% output
theta_str = num2str(theta);
rho_str = num2str(rho_d);
audiowrite(['bmt_a', num2str(theta_range(1)), '-', num2str(theta_range(2)), '_r', rho_str, '.wav'], y_out, Fs, 'BitsPerSample', 16);

