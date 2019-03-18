clear all;

% parameters
azi_density = 181;
rho_density = 100;
freq_density = 512;
rho = 1.15.^(ones(1, rho_density)+(linspace(0, 249, rho_density)/10));
theta = linspace(0, pi, azi_density);
n = 128;

% read in STF FIRs and generate SPT filters
load(['B_a', num2str(length(theta)), '_r', num2str(length(rho)), '_n', num2str(n)]);
for j=1:length(theta)
    for k=1:length(rho)
        H_FIR(j, k) = dfilt.dffir(B(j, k, :));
    end
end

% determine truncation value kk per azimuth/range point
for j=1:length(theta)
    for k=1:length(rho)
        [~, D, ~] = bmt_eigs(H_FIR(j, k));
        D = abs(diag(D));
        KK(j, k) = find((D < D(1)*0.01), 1);
        DD(:, j, k) = D;
    end
end

% reduced-order IIR modelling
for j=1:length(theta)
    for k=1:length(rho)
        H_IIR(j, k) = bmt_fun(H_FIR(j, k), KK(j, k));
    end
end

save(['IIR_a', num2str(length(theta)), '_r', num2str(length(rho)), '_n', num2str(n)], 'H_IIR');

%%

% ------------
% QC ---------
% ------------

% % KK across azimuth/range
% figure;
% surf(rho, theta*180/pi, KK);
% colormap hsv;
% xlabel('Range \rho');
% ylabel('Azimuth \theta');
% title('K_{i} variation within frontal transverse plane')
% set(gca, 'FontSize', 16);

% % KK in near and far-field
% figure;
% hold on;
% plot(theta*180/pi, KK(:,1), 'k--');
% plot(theta*180/pi, KK(:,end), 'k');
% xlabel('Azimuth \theta');
% ylabel('K_{i}');
% legend('near-field \rho = 1.25', 'far-field \rho = 37.3314');

% % HSV dropoff plots (across azimuths)
% figure;

% % normalise DD
% for i=1:azi_density
%     for j=1:rho_density
%         DD(:,i,j) = DD(:,i,j)/max(DD(:,i,j));
%     end
% end

% % near-field
% plot(DD(:, 1:20:end, 1))
% xlabel('index')
% xlim([1 10])
% ylabel('HSV')
% title('Hankel Singular Values of STF FIR models in frontal plane in the near-field (\rho = 1.15)')
% legend(sprintf('0%c', char(176)),sprintf('20%c', char(176)),sprintf('40%c', char(176)),sprintf('60%c', char(176)),sprintf('80%c', char(176)),sprintf('100%c', char(176)),sprintf('120%c', char(176)),sprintf('140%c', char(176)),sprintf('160%c', char(176)),sprintf('180%c', char(176)));
% set(gca, 'FontSize', 16);

% % far-field
% plot(DD(:, 1:20:end, end))
% xlabel('index')
% xlim([1 10])
% ylabel('HSV')
% title('Hankel Singular Values of STF FIR models in frontal plane in the far-field (\rho = 37.3314)')
% legend(sprintf('0%c', char(176)),sprintf('20%c', char(176)),sprintf('40%c', char(176)),sprintf('60%c', char(176)),sprintf('80%c', char(176)),sprintf('100%c', char(176)),sprintf('120%c', char(176)),sprintf('140%c', char(176)),sprintf('160%c', char(176)),sprintf('180%c', char(176)));
% set(gca, 'FontSize', 16);

% compare IIR response with STF response
freq_density = 512;
f = linspace(0.5,24000,freq_density);
load(['STF_H_f', num2str(length(f)), '_a', num2str(length(theta)), '_r', num2str(length(rho))]);
IIR_hh = zeros(length(f), length(theta), length(rho));
for j=1:length(theta)
    for k=1:length(rho)
        [IIR_hh(:, j, k), w] = freqz(H_IIR(j, k), f, 2*max(f));
    end
end

% % near-field
% figure;
% for j=1:20:length(theta)
%     semilogx(f, 20*log10(abs(H(:, j, 1))), 'k')
%     hold on
%     semilogx(f, 20*log10((abs(IIR_hh(:, j, 1)))), '-.r');
%     hold on
% end
% tit = strcat('STF and DBMT near-field magnitude response whilst varying \theta: (\rho = 1.15)');
% title(tit)
% xlim([100 20000])
% ylim([-40 30])
% xlabel('frequency (Hz)')
% ylabel('Magnitude (dB)')
% legend('STF', 'BMT')
% set(gca, 'FontSize', 16);
% grid on;

% % far-field
% figure;
% for j=1:20:length(theta)
%     semilogx(f, 20*log10(abs(H(:, j, end))), 'k')
%     hold on
%     semilogx(f, 20*log10((abs(IIR_hh(:, j, end)))), '-.r');
%     hold on
% end
% tit = strcat('STF and DBMT far-field magnitude response whilst varying \theta: (\rho = 37.3314)');
% title(tit)
% xlim([100 20000])
% ylim([-40 30])
% xlabel('frequency (Hz)')
% ylabel('Magnitude (dB)')
% legend('STF', 'BMT')
% set(gca, 'FontSize', 16);
% grid on;

% LSD plot
figure;
band = [50 22000];
[~, fmin] = min(abs(f-band(1)));
[~, fmax] = min(abs(f-band(2)));
f_band = f(fmin:fmax);
H_band = H(fmin:fmax, :, :);
H_IIR_band = IIR_hh(fmin:fmax, :, :);
N = length(f(fmin:fmax));
LSD = zeros(length(theta), length(rho));
for k=1:length(rho)
    for j=1:length(theta)
        LSD(j, k) = sqrt((1/N)*(sum(20*log10(abs(H_band(:, j, k))./abs(H_IIR_band(:, j, k)))).^2));
    end
end
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
title('LSD of DBMT')
set(gca, 'FontSize', 16);
