clearvars;close all;clc;

% General controls
optionPlot = 1;
near_field = 0;

%-------------------------------------------------------------------------%
% General params
%-------------------------------------------------------------------------%
n            = 16;      % Length of generated IRs
Fs           = 48000;   % Sample rate (full bandwidth)
a_h          = 0.0875;  % Head radius
c            = 343;     % Speed of sound
azi_density  = 91;      % Azimuthal sampling density
rho_density  = 10;      % Distance sampling density
freq_density = 100;      % Frequency sampling density

%-------------------------------------------------------------------------%
% Spatial/frequency vectors
%-------------------------------------------------------------------------%
rho     = 1.15.^(ones(1, rho_density)+(linspace(0, 249, rho_density)/10));
f       = linspace(0, Fs/2, freq_density);
theta   = linspace(0, pi, azi_density);

%-------------------------------------------------------------------------%
% Load STF (returns "H" matrix that contains STF's complex freq responses)
%-------------------------------------------------------------------------%
load(['STF_H_f', num2str(length(f)), '_a', num2str(length(theta)), '_r', num2str(length(rho))]);

%-------------------------------------------------------------------------%
% Derived parameters
%-------------------------------------------------------------------------%
r   = rho*a_h;
mu  = 2*pi*f*a_h/c;
f_n = (f'/max(f));                              % Normalised freq vec for "fir2" (0:1 over 0:Fs/2)

%-------------------------------------------------------------------------%
% STF --> FIR
%-------------------------------------------------------------------------%
B = zeros(length(theta), length(rho), n+1);     % Empty matrix for FIR coeffs
B_CS = B;
B_RF = B;

% Run through all azimiths and distances and generate filter approximations
for j=1:length(theta)
    for k=1:length(rho)
        % Use Matlab's "fir2" to fit the requested MAGNITUDE response
        %   --> abs(H) is sampled at the normalised freqs "f_n" (0:Fs/2)
        %   --> Return "n" feed-forward coeffs 
        b = fir2(n, f_n, abs(H(:, j, k)));
                        
        % Try to convert returned FIR taps to minimum phase
        b_CS = circshift(b, (n/2));    % simple shift (why doesn't this work?)
        b_RF = lp_fir_2_mp_fir(b);     % root-finding
                
        B(j, k, :)      = b;
        B_CS(j, k, :)   = b_CS;
        B_RF(j, k, :)   = b_RF;
        
        H_FIR(j, k)     = dfilt.dffir(b);
        H_FIR_CS(j, k)  = dfilt.dffir(b_CS);
        H_FIR_RF(j, k)  = dfilt.dffir(b_RF);
        
        %clear b
    end
end
%save(['B_a', num2str(length(theta)), '_r', num2str(length(rho)), '_n', num2str(n)], 'B');


% Compare returned frequency responses
thisThetaInd = 5;
thisRhoInd   = 5;
fig3 = figure(3);
hold on
grid on
plot(20*log10(abs(fft(H_FIR(thisThetaInd,thisRhoInd).Numerator))))
plot(20*log10(abs(fft(H_FIR_CS(thisThetaInd,thisRhoInd).Numerator))))
plot(20*log10(abs(fft(H_FIR_RF(thisThetaInd,thisRhoInd).Numerator))))





%-------------------------------------------------------------------------%
% Chris's plotting stuff
%-------------------------------------------------------------------------%

% % freq. response of FIR bank across range/azimuth
% FIR_hh = zeros(length(f), length(theta), length(rho));
% for j=1:length(theta)
%     for k=1:length(rho)
%         [FIR_hh(:, j, k), w] = freqz(H_FIR(j, k), f, 2*max(f));
%     end
% end
% 
% % compare FIR response with STF response
% if optionPlot
%     figure;
%     if near_field
%         % near-field
%         for j=1:10:length(theta)
%             %semilogx(mu, 20*log10(abs(H(:, j, 1))), 'k')
%             plot(mu, 20*log10(abs(H(:, j, 1))), 'k')
%             hold on
%             %semilogx(mu, 20*log10((abs(FIR_hh(:, j, 1)))), '-.r');
%             plot(mu, 20*log10((abs(FIR_hh(:, j, 1)))), '-.r');
%             %hold on
%         end
%         tit = strcat('STF and FIR magnitude response whilst varying AOI (\theta): \rho = 1.25');
%         title(tit)
%         xlim([0.1 35])
%         ylim([-40 30])
%         xlabel('\mu = 2\pifa/c (normalized frequency)')
%         ylabel('Magnitude (dB)')
%         legend('STF', 'BMT')
%         set(gca, 'FontSize', 16);
%     
%     else        
%         % far-field
%         for j=1:10:length(theta)
%             %semilogx(mu, 20*log10(abs(H(:, j, end))), 'k')
%             plot(mu, 20*log10(abs(H(:, j, end))), 'k')
%             hold on
%             %semilogx(mu, 20*log10((abs(FIR_hh(:, j, end)))), '-.r');
%             plot(mu, 20*log10((abs(FIR_hh(:, j, end)))), '-.r');
%             hold on
%         end
%         tit = strcat('STF and FIR magnitude response whilst varying AOI (\theta): \rho = 37.3314');
%         title(tit)
%         %xlim([0.1 35])
%         %xlim([10^(-4) 10^(2)])
%         ylim([-25 10])
%         xlabel('\mu = 2\pifa/c (normalized frequency)')
%         ylabel('Magnitude (dB)')
%         legend('STF', 'BMT')
%         set(gca, 'FontSize', 16);
%     end
% end