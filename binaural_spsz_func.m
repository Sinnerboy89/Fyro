function y = binaural_spsz_func(x,gamma,Fs,c,a_h)
%BINAURAL SPSZ Filters input with the single-pole single-zero structural
%binaural model based on Brown & Duda's 1998 paper [cite]. Note input
%signal needs to be in column vector format.

%% Internal parameters and initializations

gamma_ears = [-pi/2, pi/2]; % azimuths of left and right ear
e2e = a_h*2; % distance between ears, used for ITDs
alpha_min = 0.1; % values chosen by B & D as reasonable starting points before personalization
theta_min = 150*2*pi/360; % values chosen by B & D as reasonable starting points before personalization
L = length(x);

w0 = c/a_h;

y = zeros(L,2);

%% AOIs and ITD calculations 
 
if sign(gamma)==-1 % i.e. source is in front-left quarter space
    theta = [-gamma_ears(1)-abs(gamma),gamma_ears(2)+abs(gamma)]; % left ear is ipsilateral, right ear is contralateral
    Td = [-(e2e/(2*c))*cos(theta(1))+(e2e/(2*c)),(e2e/(2*c))*(abs(theta(2))-(pi/2))+(e2e/(2*c))]; % time delays between centre and left and right ears
else % i.e. source is in front-right quarter space
    theta = [-gamma_ears(1)+abs(gamma),gamma_ears(2)-abs(gamma)]; % vice versa
    Td = [(e2e/(2*c))*(abs(theta(1))-(pi/2))+(e2e/(2*c)),-(e2e/(2*c))*cos(theta(2))+(e2e/(2*c))];
end

Td_n = Td*Fs; % delay in 'floating sample points' for left and right ear

%% Spherical frequency response approximation
%alpha represents an approximation to the frequency response of an ideal
%rigid sphere (as described by Rayleigh), and is valid for
%0<abs(gamma)<5*pi/6. Beyond this angle, response starts to rise to the
%"bright spot" at pi - but this is non-linear, and not simulated here.

alpha = [(1+alpha_min/2)+(1-alpha_min/2)*cos(theta(1)*pi/theta_min),(1+alpha_min/2)+(1-alpha_min/2)*cos(theta(2)*pi/theta_min)];

%% Head shadow IIRs

B = [w0+(alpha*Fs);w0-(alpha*Fs)];
a = [w0+Fs,w0-Fs];

HS_IIR = [dfilt.df1(B(:,1),a),dfilt.df1(B(:,2),a)];

%% Signal flow

y = [filter(HS_IIR(1),x),filter(HS_IIR(2),x)]; % ILD IIR
y_l = GHZ_fracDelay(y(:,1)', Td_n(1))'; % ITD delay
y_r = GHZ_fracDelay(y(:,2)', Td_n(2))'; % ITD delay
if length(y_l)>length(y_r)
    y_r = [y_r;zeros(length(y_l)-length(y_r),1)];
elseif length(y_r)>length(y_l)
    y_l = [y_l;zeros(length(y_r)-length(y_l),1)]; % zero-padding shorter channel
end
y = [y_l,y_r];

end

