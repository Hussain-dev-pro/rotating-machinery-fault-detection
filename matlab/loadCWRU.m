function [x, rpm, msg] = loadCWRU(num, datadir, fs, rpmNom, checkRate)
%LOADCWRU  Drive-end signal and speed from a CWRU .mat file.
%   Picks variable Xnnn_DE_time (or the first *_DE_time), reads *RPM if present,
%   otherwise uses rpmNom. If checkRate is true, tests whether the file was really
%   recorded at 4x the assumed rate (shaft harmonics would appear at fr/4) and if so
%   low-pass filters and decimates it to fs.
S  = load(fullfile(datadir, sprintf('%d.mat', num)));  nm = fieldnames(S);
de = sprintf('X%03d_DE_time', num);
if ~isfield(S, de), de = nm{find(endsWith(nm,'_DE_time'),1)}; end
x = double(S.(de)(:));
r = nm(endsWith(nm,'RPM'));
if ~isempty(r) && ~isempty(S.(r{1})), rpm = double(S.(r{1})(1)); else, rpm = rpmNom; end
msg = '';
if nargin > 4 && checkRate
    fr = rpm/60;  n = min(numel(x), 8*fs);  y = x(1:n) - mean(x(1:n));
    A = abs(fft(y));  f = (0:n-1)'*fs/n;
    pk = @(f0) max(A(f > 0.97*f0 & f < 1.03*f0)) / median(A(f > 0.5*f0 & f < 1.5*f0));
    ratio = (pk(fr/4) + pk(fr/2) + pk(3*fr/4)) / (pk(fr) + pk(2*fr) + pk(3*fr));
    if ratio > 1.5
        n = numel(x);  X = fft(x);  f = (0:n-1)'*4*fs/n;
        X(f > 0.42*fs & f < 4*fs - 0.42*fs) = 0;  x = real(ifft(X));  x = x(1:4:end);
        msg = sprintf('shaft-line ratio %.2f -> recorded at %d kHz, decimated to %d kHz', ratio, 4*fs/1000, fs/1000);
    else
        msg = sprintf('shaft-line ratio %.2f -> recorded at %d kHz', ratio, fs/1000);
    end
end
end
