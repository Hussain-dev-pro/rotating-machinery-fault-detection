%% RMF_Explore.m - design study: which indicator / classifier generalises to unseen loads?
% Run after Phases 1-2 of RMF_Main (needs Fm, Xs in the workspace).
fs = 12000;  y = Fm(:,2);  tr = Fm(:,4) <= 1;  te = ~tr;  N = size(Xs,1);
RA = zeros(N,3); RB = RA; RC = RA; RD = RA;
for j = 1:N
    [fe,Ae,~,~,~,Ar] = envSpectrum(Xs(j,:)', fs);  Bf = bearingFreqs(Fm(j,5));
    RA(j,:) = lpInd(fe,Ae,Bf,[1 2 4]);  RB(j,:) = lpInd(fe,Ae,Bf,[2 4]);
    RC(j,:) = lpInd(fe,Ae,Bf,[1 2]);    RD(j,:) = lpInd(fe,Ar,Bf,[1 2 4]);
end
V = {Fm(:,9:11),'v1 global floor'; RA,'local, ball 1,2,4'; RB,'local, ball 2,4'; RC,'local, ball 1,2'; RD,'local, RAW spectrum'};
D = {zeros(N,0),'none'; Fm(:,7),'kurtosis'; Fm(:,6),'RMS'; [Fm(:,7) Fm(:,6)],'kurt+RMS'};
fprintf('%-20s %-10s  train  test\n','indicator','detector');
for a = 1:size(V,1), for b = 1:size(D,1)
    [P, ac] = train2(V{a,1}(tr,:), D{b,1}(tr,:), y(tr));
    bt = balAcc(y(te), classify2(V{a,1}(te,:), D{b,1}(te,:), P));
    fprintf('%-20s %-10s  %5.1f  %5.1f\n', V{a,2}, D{b,2}, 100*ac, 100*bt);
end, end
for c = 1:4, k = y==c; fprintf('%-11s RA %5.1f %5.1f %5.1f  K %.1f RMS %.3f\n', cls{c}, mean(RA(k,:)), mean(Fm(k,7)), mean(Fm(k,6))); end
