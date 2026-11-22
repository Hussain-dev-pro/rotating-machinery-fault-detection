function R = lpInd(f, A, Bf, hb)
%LPIND  Local-prominence fault indicators (dB) from an (envelope) spectrum.
%   For each harmonic h of a defect frequency F the peak inside +/-1 % of h*F
%   is compared with the MEDIAN of the spectrum in the surrounding +/-15 %
%   (the local noise floor). Averaging these dB values over the harmonics gives
%   one indicator per fault type: [outer inner ball].
%   Races use harmonics 1,2,3. Ball uses hb x BSF (default [1 2 4]); 3 x BSF is
%   skipped because it sits on the 7th shaft harmonic (7.07 x fr).
if nargin < 4, hb = [1 2 4]; end
H = {[1 2 3]*Bf.BPFO, [1 2 3]*Bf.BPFI, hb*Bf.BSF};  R = zeros(1,3);
for j = 1:3
    v = zeros(1, numel(H{j}));
    for q = 1:numel(H{j})
        F0 = H{j}(q);
        pk = max(A(f > 0.99*F0 & f < 1.01*F0));
        lo = A((f > 0.85*F0 & f < 0.97*F0) | (f > 1.03*F0 & f < 1.15*F0));
        v(q) = 20*log10(pk/median(lo));
    end
    R(j) = mean(v);
end
end
