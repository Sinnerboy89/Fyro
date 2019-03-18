clearvars;close all;clc;

% General controls
optionPlot = 1;
near_field = 0;

%-------------------------------------------------------------------------%
% General params
%-------------------------------------------------------------------------%
n            = 32;      % Length of generated IRs
Fs           = 48000;   % Sample rate (full bandwidth)
a_h          = 0.0875;  % Head radius
c            = 343;     % Speed of sound
azi_density  = 91;      % Azimuthal sampling density
rho_density  = 10;      % Distance sampling density
freq_density = 128;      % Frequency sampling density

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
f_n = (f'/max(f));  % Normalised freq vec for "fir2" (0:1 over 0:Fs/2)

%-------------------------------------------------------------------------%
% STF --> FIR
%-------------------------------------------------------------------------%
B = zeros(length(theta), length(rho), n+1); % Empty matrix for FIR coeffs
B_CS = B;
B_RF = B;

% Run through all azimiths and distances and generate filter approximations
colVals = get(groot,'DefaultAxesColorOrder');
figfont = 14;

fig3 = figure(3);
fig3.Position = [39 127 1361 678];
ax(1) = subplot(3,2,1);
grid(ax(1),'on');
pp1 = animatedline(ax(1),'Color',colVals(1,:),'LineWidth',2);
pp2 = animatedline(ax(1),'Color',colVals(2,:),'LineWidth',2); 
pp3 = animatedline(ax(1),'Color',colVals(3,:),'LineWidth',2);
ax(1).Title.String = ['FIR Impulse Response'];
ax(1).XLabel.String = ['Time (index)'];

ax(2) = subplot(3,2,3);
grid(ax(2),'on');
pp4 = animatedline(ax(2),'Color',colVals(1,:),'LineWidth',2);
pp5 = animatedline(ax(2),'Color',colVals(2,:),'LineWidth',2); 
pp6 = animatedline(ax(2),'Color',colVals(3,:),'LineWidth',2);
ax(2).Title.String = ['Magnitude response (curves offset)'];
ax(2).XLabel.String = ['Frequency (bin)'];

ax(3) = subplot(3,2,5);
grid(ax(3),'on');
pp7 = animatedline(ax(3),'Color',colVals(1,:),'LineWidth',2);
pp8 = animatedline(ax(3),'Color',colVals(2,:),'LineWidth',2); 
pp9 = animatedline(ax(3),'Color',colVals(3,:),'LineWidth',2);
ylim(ax(3),[-50,50])
ax(3).Title.String = ['Phase response (unwrapped)'];
ax(3).XLabel.String = ['Frequency (bin)'];
ax(3).YLabel.String = ['Phase (radians)'];

% Set up the Pole-Zero axes
ax(4) = subplot(3,2,2);
axis square
ax(4).Title.String = ['STF-->FIR'];
ax(5) = subplot(3,2,4);
axis square
ax(5).Title.String = ['STF-->FIR-->Rootfind'];
ax(6) = subplot(3,2,6);
axis square
ax(6).Title.String = ['STF-->FIR-->Circshift'];

for j=1:length(theta)
%for j=1:1
    for k=1:length(rho)
    %for k=1:1
        for nShift = 1:round(n/2)+1
            % Use Matlab's "fir2" to fit the requested MAGNITUDE response
            %   --> abs(H) is sampled at the normalised freqs "f_n" (0:Fs/2)
            %   --> Return "n" feed-forward coeffs
            b = fir2(n, f_n, abs(H(:, j, k)));
            
            % Try to convert returned FIR taps to minimum phase           
            b_RF = lp_fir_2_mp_fir(b);
            b_CS = circshift(b, nShift);
            %b_CS = circshift(b, (n/2)+1);
            
            B(j, k, :)      = b;
            B_CS(j, k, :)   = b_CS;
            B_RF(j, k, :)   = b_RF;
            
            H_FIR(j, k)     = dfilt.dffir(b);            
            H_FIR_RF(j, k)  = dfilt.dffir(b_RF);
            H_FIR_CS(j, k)  = dfilt.dffir(b_CS);
            
            
            % Compare returned frequency responses
            thisThetaInd = j;
            thisRhoInd   = k;
                        
            clearpoints(pp1);clearpoints(pp2);clearpoints(pp3);
            addpoints(pp1,[1:n+1],permute(B(thisThetaInd,thisRhoInd,:),[3,2,1]));            
            addpoints(pp2,[1:n+1],permute(B_RF(thisThetaInd,thisRhoInd,:),[3,2,1]));            
            addpoints(pp3,[1:n+1],permute(B_CS(thisThetaInd,thisRhoInd,:),[3,2,1]));
            l1 = legend(ax(1),['STF-->FIR'],['STF-->FIR-->Rootfind'],['STF-->FIR-->Circshift']);
            l1.FontSize = figfont-2;
                                    
            clearpoints(pp4);clearpoints(pp5);clearpoints(pp6);
            addpoints(pp4,[1:n+1],20*log10(abs(fft(H_FIR(thisThetaInd,thisRhoInd).Numerator))));            
            addpoints(pp5,[1:n+1],2+20*log10(abs(fft(H_FIR_RF(thisThetaInd,thisRhoInd).Numerator))));            
            addpoints(pp6,[1:n+1],1+20*log10(abs(fft(H_FIR_CS(thisThetaInd,thisRhoInd).Numerator))));                                    
            
            clearpoints(pp7);clearpoints(pp8);clearpoints(pp9);
            addpoints(pp7,[1:n+1],unwrap(angle(fft(H_FIR(thisThetaInd,thisRhoInd).Numerator))));
            addpoints(pp8,[1:n+1],unwrap(angle(fft(H_FIR_RF(thisThetaInd,thisRhoInd).Numerator))));
            addpoints(pp9,[1:n+1],unwrap(angle(fft(H_FIR_CS(thisThetaInd,thisRhoInd).Numerator))));
            
            cla(ax(4));cla(ax(5));cla(ax(6));
            zplane(b,1,ax(4))
            zplane(b_RF,1,ax(5))
            zplane(b_CS,1,ax(6))                        
            
            ax(2).Title.String = {['Magnitude response (curves offset). Azimuth: ' num2str(theta(j)) '; Distance: ' num2str(rho(k))];...
                ['Number of samples of circshift: ' num2str(nShift)]};
            
            drawnow;
        end
    end
end
%save(['B_a', num2str(length(theta)), '_r', num2str(length(rho)), '_n', num2str(n)], 'B');






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