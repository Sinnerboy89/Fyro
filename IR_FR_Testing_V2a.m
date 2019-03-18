% Testing the notion of IR and FR
clearvars;close all;clc

% Number of points to use in plotting frequency response (interpolation!)
N_Points = 512*2;
% N_Points = 64;

% Make a random impulse response, then zero pad it
IR      = rand(64,1);
IR      = [IR;zeros(N_Points-length(IR),1)];

% Make another version of the IR that is circularly shifted
IR_v2   = circshift(IR,round(length(IR)/2));
%IR_v2   = circshift(IR,10);

% Make FRs via FFT
FR      = fft(IR);
FR_v2   = fft(IR_v2);

% Make FRs via freqz (these should be exactly the same as the ones via FFT)
[FR_freqz,omega]    = freqz(IR,1,N_Points,'whole');
[FR_v2_freqz,~]     = freqz(IR_v2,1,N_Points,'whole');


%-------------------------------------------------------------------------%
% Plotting
%-------------------------------------------------------------------------%
fig1 = figure(1);
subplot(3,2,1)
plot(IR)
hold on
grid on
plot(IR_v2)
legend('Zero-padded IR','Zero-padded+circularly shifted IR')
xlabel('Time (samples)')
ylabel('Amplitude')

subplot(3,2,3)
plot(omega/(2*pi),20*log10(abs(FR)))
hold on
grid on
plot(omega/(2*pi),20*log10(abs(FR_freqz)))
legend('Zero-padded IR via FFT','Zero-padded IR via freqz')
xlabel('Frequency (normalised)')
ylabel('Magnitude')

subplot(3,2,5)
plot(omega/(2*pi),20*log10(abs(FR_v2)))
hold on
grid on
plot(omega/(2*pi),20*log10(abs(FR_v2_freqz)))
legend('Zero-padded+circularly shifted IR via FFT','Zero-padded+circularly shifted IR via freqz')
xlabel('Frequency (normalised)')
ylabel('Magnitude')

subplot(3,2,4)
plot(omega/(2*pi),unwrap(angle(FR)))
hold on
grid on
plot(omega/(2*pi),unwrap(angle(FR_freqz)))
legend('Zero-padded IR via FFT','Zero-padded IR via freqz')
xlabel('Frequency (normalised)')
ylabel('Phase (radians)')

subplot(3,2,6)
plot(omega/(2*pi),unwrap(angle(FR_v2)))
hold on
grid on
plot(omega/(2*pi),unwrap(angle(FR_v2_freqz)))
legend('Zero-padded+circularly shifted IR via FFT','Zero-padded+circularly shifted IR via freqz')
xlabel('Frequency (normalised)')
ylabel('Phase (radians)')