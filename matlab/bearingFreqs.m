function F = bearingFreqs(fr)
%BEARINGFREQS  Defect frequencies (Hz) of the CWRU drive-end bearing SKF 6205-2RS JEM
%   F = bearingFreqs(fr), fr = shaft speed in Hz. Contact angle 0 (deep groove).
%   F.BPFO outer race, F.BPFI inner race, F.BSF ball spin, F.FTF cage, F.fr shaft.
%   Geometry from the CWRU Bearing Data Center (inches): 9 balls, ball diameter
%   0.3126, pitch diameter 1.537.
n = 9;  d = 0.3126;  D = 1.537;  r = d/D;
F.fr   = fr;
F.FTF  = fr/2*(1 - r);
F.BPFO = n*fr/2*(1 - r);
F.BPFI = n*fr/2*(1 + r);
F.BSF  = fr*D/(2*d)*(1 - r^2);
end
