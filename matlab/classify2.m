function yh = classify2(R, Dt, P)
%CLASSIFY2  Two-stage rule-based diagnosis.
%  R  = [Ro Ri Rb] indicators (dB), Dt = detection features (e.g. kurtosis, RMS)
%  P  = [To Ti Tb Td...] thresholds.
%  1) race fault if the larger race margin (Ro-To or Ri-Ti) is positive
%  2) otherwise ball fault if Rb > Tb or any detection feature exceeds its threshold
%     (= fault energy present but no race signature)
%  3) otherwise healthy.
n = size(R,1);  yh = ones(n,1);
[m, c] = max(R(:,1:2) - P(1:2), [], 2);
race = m > 0;
yh(race & c == 1) = 4;  yh(race & c == 2) = 2;
det = R(:,3) > P(3);
if size(Dt,2) > 0, det = det | any(Dt > P(4:end), 2); end
yh(~race & det) = 3;
end
