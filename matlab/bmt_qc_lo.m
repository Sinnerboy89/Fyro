clear all;

% parameters
azi_density = 91;
rho_density = 10;
freq_density = 128;
a_h = 0.0875;
c = 343;
rho = 1.15.^(ones(1, rho_density)+(linspace(0, 249, rho_density)/10));
f = linspace(0.1, 24000, freq_density);
theta = linspace(0, pi, azi_density);
n = freq_density*2;
near_field = 1;

% intermediate variables
Fs = 2*max(f);
r = rho*a_h;
mu = 2*pi*f*a_h/c;

% read in STF
load(['STF_H_f', num2str(length(f)), '_a', num2str(length(theta)), '_r', num2str(length(rho))]);

% read in STF FIRs and generate SPT filters
load(['B_a', num2str(length(theta)), '_r', num2str(length(rho)), '_n', num2str(n)]);
for j=1:length(theta)
    for k=1:length(rho)
        H_FIR(j, k) = dfilt.dffir(B(j, k, :));
    end
end

% freq. response of FIR bank across range/azimuth
FIR_hh = zeros(length(f), length(theta), length(rho));
for j=1:length(theta)
    for k=1:length(rho)
        [FIR_hh(:, j, k), w] = freqz(H_FIR(j, k), f, 2*max(f));
    end
end

% compare FIR response with STF response
figure;
if near_field
    % near-field
    for j=1:10:length(theta)
        semilogx(mu, 20*log10(abs(H(:, j, 1))), 'k')
        hold on
        semilogx(mu, 20*log10((abs(FIR_hh(:, j, 1)))), '-.r');
        hold on
    end
    tit = strcat('STF and FIR magnitude response whilst varying AOI (\theta): \rho = 1.25');
    title(tit)
    xlim([0.1 35])
    ylim([-40 30])
    xlabel('\mu = 2\pifa/c (normalized frequency)')
    ylabel('Magnitude (dB)')
    legend('STF', 'BMT')
    set(gca, 'FontSize', 16);
else
    % far-field
    for j=1:10:length(theta)
        semilogx(mu, 20*log10(abs(H(:, j, end))), 'k')
        hold on
        semilogx(mu, 20*log10((abs(FIR_hh(:, j, end)))), '-.r');
        hold on
    end
    tit = strcat('STF and FIR magnitude response whilst varying AOI (\theta): \rho = 37.3314');
    title(tit)
    xlim([0.1 35])
    ylim([-25 10])
    xlabel('\mu = 2\pifa/c (normalized frequency)')
    ylabel('Magnitude (dB)')
    legend('STF', 'BMT')
    set(gca, 'FontSize', 16);
end

% determine truncation value kk per azimuth/range point
for j=1:length(theta)
    for k=1:length(rho)
        [~, D, ~] = bmt_eigs(H_FIR(j, k));
        D = abs(diag(D));
        KK(j, k) = find((D < 0.3), 1);
        DD(:, j, k) = D;
    end
end

% KK across azimuth/range
figure;
surf(rho, theta*180/pi, KK);
xlabel('Range \rho');
ylabel('Azimuth \theta');
title('K_{i} variation within frontal transverse plane')

% KK in near and far-field
figure;
hold on;
plot(theta*180/pi, KK(:,1), 'k--');
plot(theta*180/pi, KK(:,end), 'k');
xlabel('Azimuth \theta');
ylabel('K_{i}');
legend('near-field \rho = 1.25', 'far-field \rho = 37.3314');

% HSV dropoff plots (across azimuths)
figure;
if near_field
    plot(DD(:, :, 1))
    xlabel('index')
%     xlim([1 10])
    ylabel('HSV')
    title('Hankel Singular Values of STF FIR models in frontal plane in the near-field (\rho = 1.25)')
else
    plot(DD(:, :, end))
    xlabel('index')
    xlim([1 10])
    ylabel('HSV')
    title('Hankel Singular Values of STF FIR models in frontal plane in the far-field (\rho = 37.3314)')
end

% reduced-order IIR modelling
for j=1:length(theta)
    for k=1:length(rho)
        H_IIR(j, k) = bmt_fun(H_FIR(j, k), KK(j, k));
    end
end

% freq. response of IIR bank across range/azimuth
IIR_hh = zeros(length(f), length(theta), length(rho));
for j=1:length(theta)
    for k=1:length(rho)
        [IIR_hh(:, j, k), w] = freqz(H_IIR(j, k), f, 2*max(f));
    end
end

% compare IIR response with STF response
figure;
if near_field
    % near-field
    for j=1:10:length(theta)
        semilogx(mu, 20*log10(abs(H(:, j, 1))), 'k')
        hold on
        semilogx(mu, 20*log10((abs(IIR_hh(:, j, 1)))), '-.r');
        hold on
    end
    tit = strcat('STF and IIR magnitude response whilst varying AOI (\theta): \rho = 1.25');
    title(tit)
    xlim([0.1 35])
    ylim([-40 30])
    xlabel('\mu = 2\pifa/c (normalized frequency)')
    ylabel('Magnitude (dB)')
    legend('STF', 'BMT')
    set(gca, 'FontSize', 16);
else
    % far-field
    for j=1:10:length(theta)
        semilogx(mu, 20*log10(abs(H(:, j, end))), 'k')
        hold on
        semilogx(mu, 20*log10((abs(IIR_hh(:, j, end)))), '-.r');
        hold on
    end
    tit = strcat('STF and IIR magnitude response whilst varying AOI (\theta): \rho = 37.3314');
    title(tit)
    xlim([0.1 35])
    ylim([-25 10])
    xlabel('\mu = 2\pifa/c (normalized frequency)')
    ylabel('Magnitude (dB)')
    legend('STF', 'BMT')
    set(gca, 'FontSize', 16);
end

% LSD plot
band = [10 20000];
[~, fmin] = min(abs(f-band(1)));
[~, fmax] = min(abs(f-band(2)));
f_band = f(fmin:fmax);
mu_band = mu(fmin:fmax);
H_band = H(fmin:fmax, :, :);
H_BMT_band = IIR_hh(fmin:fmax, :, :);
N = length(f(fmin:fmax));
LSD = zeros(length(theta), length(rho));
for k=1:length(rho)
    for j=1:length(theta)
        LSD(j, k) = sqrt((1/N)*(sum(20*log10(abs(H_band(:, j, k))./abs(H_BMT_band(:, j, k)))).^2));
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
caxis([0 10])
colormap jet
shading interp
view(2)
title('Performance of BMT approximation to STF')
