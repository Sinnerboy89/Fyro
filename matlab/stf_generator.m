clear all;

azi_density = 181;
rho_density = 100;
freq_density = 512;

a_h = 0.0875;
rho = 1.15.^(ones(1,rho_density)+(linspace(0,249,rho_density)/10));
r = rho*a_h;
threshold = 0.001;
f = linspace(0.5,24000,freq_density);
theta = linspace(0,pi,azi_density);
% theta = [0, linspace(10*pi/180,170*pi/180,azi_density)]; # for QC

H = zeros(length(f),length(theta),length(rho));
for j=1:length(theta)
    for k=1:length(rho)
        for i=1:length(f)
            H(i,j,k) = sphere(a_h,r(k),theta(j),f(i),threshold);
        end
    end
end

save(['STF_H_f', num2str(freq_density), '_a', num2str(azi_density), '_r', num2str(rho_density)], 'H');

%%

% ------------
% QC ---------
% ------------

% c = 343;
% mu = 2*pi*f*a_h/c;

% compare responses at ipsilateral and near-contralateral locations
% figure;
% for i=1:length(rho)
%     semilogx(mu, 20*log10(abs(H(:, 1, i))), 'k', 'LineWidth', 1.25)
%     hold on;
% end
% for i=1:length(rho)
%     semilogx(mu, 20*log10(abs(H(:, end-5, i))), 'r-.', 'LineWidth', 1.25)
%     hold on;
% end
% xlim([0.1 38])
% ylim([-40 25])
% xlabel('\mu = 2\pifa/c')
% ylabel('Magnitude (dB)')
% legend({'\theta = 0', '\theta = 170'}, 'Location', 'southwest');
% set(gca, 'FontSize', 16);

% % look at response as azimuth varies (near or far)
% figure;

% % near-field
% for j=1:length(theta)
%     semilogx(f, 20*log10(abs(H(:, j, 1))))
%     hold on
% end
% tit = strcat('STF near-field magnitude response whilst varying \theta: (\rho = 1.15, a = 8.75cm)');
% title(tit)
% xlim([100 20000])
% ylim([-40 25])
% xlabel('frequency (Hz)')
% ylabel('Magnitude (dB)')
% legend(sprintf('0%c', char(176)),sprintf('10%c', char(176)),sprintf('30%c', char(176)),sprintf('50%c', char(176)),sprintf('70%c', char(176)),sprintf('90%c', char(176)),sprintf('110%c', char(176)),sprintf('130%c', char(176)),sprintf('150%c', char(176)),sprintf('170%c', char(176)));
% set(gca, 'FontSize', 16);
% grid on;

% % far-field
% for j=1:length(theta)
%     semilogx(f, 20*log10(abs(H(:, j, end))))
%     hold on
% end
% tit = strcat('STF far-field magnitude response whilst varying \theta: (\rho = 37.3314, a = 8.75cm)');
% title(tit)
% xlim([100 20000])
% ylim([-40 25])
% xlabel('frequency (Hz)')
% ylabel('Magnitude (dB)')
% legend(sprintf('0%c', char(176)),sprintf('10%c', char(176)),sprintf('30%c', char(176)),sprintf('50%c', char(176)),sprintf('70%c', char(176)),sprintf('90%c', char(176)),sprintf('110%c', char(176)),sprintf('130%c', char(176)),sprintf('150%c', char(176)),sprintf('170%c', char(176)));
% set(gca, 'FontSize', 16);
% grid on;