clear all;

% parameters
azi_density = 181;
rho_density = 100;
freq_density = 512;
% n = freq_density*2;
n = round(freq_density/4);

% intermediates
rho = 1.15.^(ones(1, rho_density)+(linspace(0, 249, rho_density)/10));
f = linspace(0.1, 24000, freq_density);
theta = linspace(0, pi, azi_density);
a_h = 0.0875;
Fs = 2*max(f);
c = 343;
r = rho*a_h;
mu = 2*pi*f*a_h/c;

% read in STF
load(['STF_H_f', num2str(length(f)), '_a', num2str(length(theta)), '_r', num2str(length(rho))]);

% STF --> FIR
B = zeros(length(theta), length(rho), n+1);
f(1) = 0;
f_n = (f'/max(f));
for j=1:length(theta)
    for k=1:length(rho)
        b = fir2(n, f_n, abs(H(:, j, k)));
        b = lp_fir_2_mp_fir(b); % root-finding
%         b = circshift(b, (n/2)+1); % simple circshift
        B(j, k, :) = b;
        H_FIR_filt(j, k) = dfilt.dffir(b);
        clear b
    end
end

save(['B_a', num2str(length(theta)), '_r', num2str(length(rho)), '_n', num2str(n)], 'B');

%%

% ------------
% QC ---------
% ------------

% % freq. response of FIR bank across range/azimuth
% FIR_hh = zeros(length(f), length(theta), length(rho));
% for j=1:length(theta)
%     for k=1:length(rho)
%         %         [FIR_hh(:, j, k), w] = freqz(H_FIR(j, k), f, 2*max(f));
%                 [~, w] = freqz(H_FIR(j, k), f, 2*max(f));
%                 twosided = fft(H_FIR(j, k).Numerator);
%                 FIR_hh(:, j, k) = twosided(1:freq_density);
%     end
% end
% 
% % compare FIR response with STF response
% if plot
%     figure;
%     if near_field
%         % near-field
%         for j=1:10:length(theta)
%             semilogx(mu, 20*log10(abs(H(:, j, 1))), 'k')
%             hold on
%             semilogx(mu, 20*log10((abs(FIR_hh(:, j, 1)))), '-.r');
%             hold on
%         end
%         tit = strcat('STF and FIR magnitude response whilst varying AOI (\theta): \rho = 1.25');
%         title(tit)
%         xlim([0.1 35])
%         ylim([-40 30])
%         xlabel('\mu = 2\pifa/c (normalized frequency)')
%         ylabel('Magnitude (dB)')
%         legend('STF', 'BMT')
%         set(gca, 'FontSize', 16);
%     else
%         % far-field
%         for j=1:10:length(theta)
%             semilogx(mu, 20*log10(abs(H(:, j, end))), 'k')
%             hold on
%             semilogx(mu, 20*log10((abs(FIR_hh(:, j, end)))), '-.r');
%             hold on
%         end
%         tit = strcat('STF and FIR magnitude response whilst varying AOI (\theta): \rho = 37.3314');
%         title(tit)
%         xlim([0.1 35])
%         ylim([-25 10])
%         xlabel('\mu = 2\pifa/c (normalized frequency)')
%         ylabel('Magnitude (dB)')
%         legend('STF', 'BMT')
%         set(gca, 'FontSize', 16);
%     end
% end

% LSD plot
for j=1:length(theta)
    for k=1:length(rho)
        [H_FIR(:, j, k), ~] = freqz(H_FIR_filt(j, k), f, 2*max(f));
    end
end
band = [50 22000];
[~, fmin] = min(abs(f-band(1)));
[~, fmax] = min(abs(f-band(2)));
f_band = f(fmin:fmax);
H_band = H(fmin:fmax, :, :);
H_FIR_band = H_FIR(fmin:fmax, :, :);
N = length(f(fmin:fmax));
LSD = zeros(length(theta), length(rho));
for k=1:length(rho)
    for j=1:length(theta)
        LSD(j, k) = sqrt((1/N)*(sum(20*log10(abs(H_band(:, j, k))./abs(H_FIR_band(:, j, k)))).^2));
    end
end
figure
surf(rho, theta*180/pi, LSD)
set(gca,  'XScale', 'log')
xlabel('range \rho')
xlim([1.25 35])
NumTicks = 6;
L = get(gca, 'XLim');
set(gca, 'XTick', linspace(L(1), L(2), NumTicks))
ylabel('angle of incidence \theta')
ylabel(colorbar, 'log-spectral distortion (dB)')
caxis([0 1])
colormap jet
shading interp
view(2)
title('LSD of STF after conversion to FIR (M = FD/2)')