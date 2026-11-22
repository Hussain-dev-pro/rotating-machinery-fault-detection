function [fe,Ae,band,kband,fr_,Ar] = envSpectrum(x, fs)
%ENVSPECTRUM  Envelope spectrum with automatic band selection (no toolboxes).
%   1. Candidate bands (500, 1000, 2000 Hz wide) across the spectrum.
%   2. Pick the band whose band-passed signal has the highest kurtosis
%      (impulsive bearing impacts excite a structural resonance).
%   3. Analytic signal by FFT (Hilbert transform), envelope = |z|.
%   4. Hann-windowed amplitude spectrum of the envelope.
%   Also returns the raw amplitude spectrum (fr_, Ar) for comparison.
x = x(:) - mean(x);  N = numel(x);
X = fft(x);  f = (0:N-1)'*fs/N;
best = -inf;  band = [0 fs/2];
for bw = [500 1000 2000]
    for fc = bw/2+250 : bw/4 : fs/2-bw/2-250
        m = f >= fc-bw/2 & f <= fc+bw/2;               % positive frequencies only
        z = ifft(2*X.*m);                              % analytic band-passed signal
        k = kurt(real(z));
        if k > best, best = k; band = [fc-bw/2 fc+bw/2]; end
    end
end
kband = best;
m = f >= band(1) & f <= band(2);
e = abs(ifft(2*X.*m));  e = e - mean(e);
w = hannw(N);
E = abs(fft(e.*w))*2/sum(w);  K = floor(N/2);
fe = f(1:K);  Ae = E(1:K);
R = abs(fft(x.*w))*2/sum(w);
fr_ = fe;  Ar = R(1:K);
end

function k = kurt(x)
    x = x - mean(x);  k = mean(x.^4)/mean(x.^2)^2;
end

function w = hannw(N)
    w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
end
