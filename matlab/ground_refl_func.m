function [ refl_p,r_r ] = ground_refl_func(f,k,d_path)
%GROUND_REFL_FUNC Summary of this function goes here
%   Detailed explanation goes here

sigma = 50; % corresponds to grass/pine forest, taken from unknown reference

h_source = 1.65;                   % height of source in m
h_receiver = 1.65;                 % height of receiver in m

% constants
Z_air = 415;                    % characteristic impedance of air in N - sec/m^3 at 20 C

% based on input

% pathlength for reflected wave
r_r = sqrt(d_path^2 + (h_source + h_receiver)^2);

angle = acos(d_path / r_r);

Z_ground_real = 1 + 9.08 .* (f ./ sigma).^(-0.75);
Z_ground_imaginary = 11.9 .* (f ./ sigma).^(-0.73);
Z_ground = Z_ground_real + 1i .* Z_ground_imaginary;

% plane reflection coefficient
R_p = (sin(angle) - (Z_air ./ Z_ground)) ./ (sin(angle) + (Z_air ./ Z_ground));

% second term that describes reflected wave pressure
refl_p =  R_p.* exp(-1i*k*r_r);     
   

























end

