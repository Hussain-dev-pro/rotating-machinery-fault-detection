%% RMF_PhoneTest.m  -  check a real machine at home with a phone accelerometer
% Put the phone flat on a running machine (table fan base, washing machine,
% exercise bike...). Record 20 s with phyphox "Acceleration (without g)" and
% export the CSV as phone_machine.csv. Enter the machine speed if known.
% The script finds the running-speed line (1x) and its harmonics, and gives
% the vibration velocity in mm/s (10 Hz - Nyquist band) as a simple health check.
clear; clc; close all;
file = 'phone_machine.csv';  axisCol = 4;  rpmGuess = [];     % e.g. 1200 for a fan
if isfile(file)
    D = readmatrix(file);  t = D(:,1);  a = D(:,axisCol);  src = 'MEASURED';
else
    fs0 = 400;  t = (0:1/fs0:20)';  f1 = 19.8;  rng(1);             % demo: slightly unbalanced fan
    a = 0.25*sin(2*pi*f1*t) + 0.06*sin(2*pi*2*f1*t+1) + 0.03*randn(size(t));  src = 'DEMO (no csv found)';
end
ok = ~isnan(t) & ~isnan(a);  t = t(ok);  a = a(ok);
fs = 1/median(diff(t));  tu = (t(1):1/fs:t(end))';  au = interp1(t,a,tu);  au = au - mean(au);
N = numel(au);  w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
A = abs(fft(au.*w))*2/sum(w);  A = A(1:floor(N/2));  f = (0:numel(A)-1)'*fs/N;
if isempty(rpmGuess), s = f > 3 & f < fs/2*0.9; else, s = f > 0.8*rpmGuess/60 & f < 1.2*rpmGuess/60; end
[~,i] = max(A.*s);  f1 = f(i);
h = zeros(1,3);
for k = 1:3, m = abs(f - k*f1) < 0.02*k*f1 + fs/N;  h(k) = max(A(m)); end
V = A./(2*pi*max(f,1e-9))*1e3;  V(f < 10) = 0;               % velocity amplitude, mm/s
vrms = sqrt(sum((V/sqrt(2)).^2));
fprintf('%s data, fs = %.0f Hz\n', src, fs);
fprintf('Running speed 1x = %.2f Hz (%.0f rpm)\n', f1, f1*60);
fprintf('Harmonic amplitudes 1x / 2x / 3x = %.3f / %.3f / %.3f m/s^2  (2x/1x = %.2f)\n', h, h(2)/h(1));
fprintf('Velocity (10 Hz - %.0f Hz): %.2f mm/s rms\n', fs/2, vrms);
if h(2)/h(1) < 0.5, fprintf('Pattern: dominant 1x -> typical of UNBALANCE\n');
else, fprintf('Pattern: strong 2x -> check MISALIGNMENT or looseness\n'); end
fig = figure('Color','w');  try, theme(fig,'light'); end %#ok<TRYNC>
plot(f, A, 'k'); hold on; for k = 1:3, xline(k*f1, 'r--', sprintf('%dx',k)); end
grid on; xlim([0 min(fs/2, 6*f1)]); xlabel('Frequency (Hz)'); ylabel('Acceleration (m/s^2)');
title(['Phone vibration spectrum: ' src]);
if ~isfolder('figures'), mkdir('figures'); end
exportgraphics(fig, fullfile('figures','fig6_phone_test.png'), 'Resolution', 150);
