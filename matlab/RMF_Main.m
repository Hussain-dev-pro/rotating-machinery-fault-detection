%% ========================================================================
%  RMF_Main.m  -  Rotating-machinery bearing fault detection on REAL data
%  ------------------------------------------------------------------------
%  Data: Case Western Reserve University (CWRU) Bearing Data Center,
%        12 kHz drive-end accelerometer, SKF 6205-2RS JEM bearing,
%        healthy + inner race + ball + outer race (@6:00) faults,
%        fault sizes 0.007/0.014/0.021 in, motor loads 0-3 hp (40 files).
%  Phases
%   1  Download data, bearing defect frequencies (theory vs CWRU table)
%   2  Segment signals (1 s) and extract features
%        time: RMS, kurtosis, crest factor
%        ENVELOPE and raw-FFT local-prominence indicators at BPFO / BPFI / BSF
%   3  Plots: raw signals, raw spectra, envelope spectra (why FFT alone fails)
%   4  Diagnosis: two-stage rule (RMS alarm + envelope signature), thresholds
%        trained on loads 0-1 hp, TESTED on unseen 2-3 hp; ablation v1 -> v2;
%        envelope method vs raw-FFT method, confusion matrices
%   5  Severity trend and robustness to added noise
%   6  Summary and save
%  Needs only base MATLAB (fft, websave). Uses bearingFreqs.m, envSpectrum.m,
%  loadCWRU.m, lpInd.m, classify2.m, train2.m, balAcc.m.
%  ========================================================================
clear; clc; close all;  rng(0);
figdir = 'figures';  if ~isfolder(figdir), mkdir(figdir); end
datadir = 'data';    if ~isfolder(datadir), mkdir(datadir); end
fs = 12000;  segN = 12000;  nSeg = 10;               % 1 s segments, 10 per file
cls = {'Healthy','Inner race','Ball','Outer race'};
rpmNom = [1797 1772 1750 1730];

%% 1) DATA AND THEORY ------------------------------------------------------
T = [97 1 0 0; 98 1 0 1; 99 1 0 2; 100 1 0 3];       % [file class size(mils) load(hp)]
IR = [105 106 107 108; 169 170 171 172; 209 210 211 212];
BA = [118 119 120 121; 185 186 187 188; 222 223 224 225];
OR = [130 131 132 133; 197 198 199 200; 234 235 236 237];
sz = [7 14 21];
for s = 1:3
    for L = 0:3
        T = [T; IR(s,L+1) 2 sz(s) L; BA(s,L+1) 3 sz(s) L; OR(s,L+1) 4 sz(s) L]; %#ok<AGROW>
    end
end
fprintf('==== PHASE 1: DATA ====\n%d CWRU files (12 kHz drive end): ', size(T,1));
opt = weboptions('Timeout', 120);  nd = 0;
for i = 1:size(T,1)
    fn = fullfile(datadir, sprintf('%d.mat', T(i,1)));
    if ~isfile(fn)
        websave(fn, sprintf('https://engineering.case.edu/sites/default/files/%d.mat', T(i,1)), opt);  nd = nd + 1;
    end
end
fprintf('%d downloaded, %d already present\n', nd, size(T,1)-nd);
B0 = bearingFreqs(1);
fprintf('Defect frequencies (x shaft speed)  computed | CWRU table\n');
fprintf('  BPFO outer race   %.4f | 3.5848\n  BPFI inner race   %.4f | 5.4152\n', B0.BPFO, B0.BPFI);
fprintf('  BSF  ball spin    %.4f | 4.7135 = 2 x BSF (ball hits both races)\n  FTF  cage         %.4f | 0.39828\n\n', B0.BSF, B0.FTF);

%% 2) SEGMENT AND EXTRACT FEATURES -----------------------------------------
fprintf('==== PHASE 2: FEATURES ====\n');
% columns: 1 file 2 class 3 size 4 load 5 fr 6 rms 7 kurt 8 crest
%          9-11 envelope Ro Ri Rb (local prominence)  12-14 same on raw FFT
%          15-16 band  17 band kurtosis  18-20 v1 indicators (global floor, kept for comparison)
Fm = [];  Xs = [];  info = cell(size(T,1),1);
for i = 1:size(T,1)
    [x, rpm, msg] = loadCWRU(T(i,1), datadir, fs, rpmNom(T(i,4)+1), T(i,2) == 1);
    if ~isempty(msg), fprintf('  %d.mat (healthy): %s\n', T(i,1), msg); end
    fr = rpm/60;  Bf = bearingFreqs(fr);
    info{i} = sprintf('%d.mat: %d samples (%.1f s), %.0f rpm', T(i,1), numel(x), numel(x)/fs, rpm);
    ns = min(nSeg, floor(numel(x)/segN));
    for s = 1:ns
        seg = x((s-1)*segN + (1:segN));
        [fe,Ae,band,kb,~,Ar] = envSpectrum(seg, fs);
        Fm = [Fm; i T(i,2:4) fr rmsN(seg) kurtosisN(seg) max(abs(seg))/rmsN(seg) ...
              lpInd(fe,Ae,Bf) lpInd(fe,Ar,Bf) band kb indicators(fe,Ae,Bf)]; %#ok<AGROW>
        Xs = [Xs; seg']; %#ok<AGROW>
    end
end
y = Fm(:,2);  load_ = Fm(:,4);
fprintf('%d segments of 1 s: healthy %d, inner %d, ball %d, outer %d\n', numel(y), sum(y==1), sum(y==2), sum(y==3), sum(y==4));
fprintf('Example: %s\n\n', info{1});
fprintf('  Class        kurtosis   RMS (g)   R_outer  R_inner  R_ball  (envelope, dB above local floor)\n');
for c = 1:4
    k = y == c;
    fprintf('  %-11s %8.1f %9.3f %8.1f %8.1f %7.1f\n', cls{c}, mean(Fm(k,7)), mean(Fm(k,6)), mean(Fm(k,9)), mean(Fm(k,10)), mean(Fm(k,11)));
end
fprintf('\n');

%% 3) WHY ENVELOPE ANALYSIS: example plots ---------------------------------
ex = [97 105 118 130];                                % 0 hp, 0.007 in faults
fig = newFig();  tl = tiledlayout(fig,4,3,'TileSpacing','compact','Padding','compact');
for c = 1:4
    i = find(T(:,1) == ex(c));  k = find(Fm(:,1) == i, 1);
    seg = Xs(k,:)';  Bf = bearingFreqs(Fm(k,5));
    [fe,Ae,band,~,~,Ar] = envSpectrum(seg, fs);
    nexttile(tl);  t = (0:segN-1)/fs;  plot(t(1:1200)*1e3, seg(1:1200), 'k');  grid on;
    ylabel(cls{c},'FontWeight','bold');  if c == 1, title('Time signal (first 100 ms)'); end
    if c == 4, xlabel('Time (ms)'); end
    nexttile(tl);  plot(fe/1e3, Ar, 'Color',[0.3 0.3 0.3]);  grid on;  xlim([0 6]);
    if c == 1, title('Raw FFT spectrum'); end
    if c == 4, xlabel('Frequency (kHz)'); end
    nexttile(tl);  plot(fe, Ae, 'b');  hold on;  grid on;  xlim([0 500]);
    yl = max(Ae(fe>5 & fe<500))*1.15;  ylim([0 yl]);
    mk = {Bf.BPFO,'r',Bf.BPFI,[0 .6 0],Bf.BSF,[0.8 0.4 0]};
    for q = 1:3, for h = 1:3, xline(h*mk{2*q-1}, ':', 'Color', mk{2*q}, 'LineWidth', 1.2); end, end
    if c == 1, title('Envelope spectrum (red BPFO, green BPFI, orange BSF)'); end
    if c == 4, xlabel('Frequency (Hz)'); end
    text(0.98,0.9,sprintf('band %.1f-%.1f kHz',band/1e3),'Units','normalized','HorizontalAlignment','right','FontSize',8);
end
exportgraphics(fig, fullfile(figdir,'fig1_signals_spectra_envelope.png'), 'Resolution', 150);

fig = newFig();
subplot(1,2,1); gscatter2(Fm(:,10), Fm(:,9), y, cls);  grid on;
xlabel('Inner-race indicator R_i (dB)'); ylabel('Outer-race indicator R_o (dB)');
title('Race faults: envelope indicators');
subplot(1,2,2); gscatter2(Fm(:,6), Fm(:,11), y, cls);  grid on;  set(gca,'XScale','log');
xlabel('Overall RMS (g, log scale)'); ylabel('Ball indicator R_b (dB)');
title('Healthy vs ball: overall level separates');
exportgraphics(fig, fullfile(figdir,'fig2_feature_space.png'), 'Resolution', 150);

%% 4) DIAGNOSIS: train on 0-1 hp, test on 2-3 hp ---------------------------
fprintf('==== PHASE 4: DIAGNOSIS (train 0-1 hp, test on unseen 2-3 hp) ====\n');
tr = load_ <= 1;  te = ~tr;  rm = Fm(:,6);  nul = zeros(numel(y),0);
% ablation: how each design step changes accuracy on the UNSEEN loads
abl = {'v1: global-floor envelope indicators', Fm(:,18:20), nul;
       'v2a: local-prominence envelope, no RMS stage', Fm(:,9:11), nul;
       'v2b: local-prominence RAW FFT + RMS stage', Fm(:,12:14), rm;
       'v2 : local-prominence ENVELOPE + RMS stage', Fm(:,9:11), rm};
fprintf('  Design step                                      train   test (balanced accuracy)\n');
abA = zeros(size(abl,1),2);
for a = 1:size(abl,1)
    [P, ac] = train2(abl{a,2}(tr,:), abl{a,3}(tr,:), y(tr));
    abA(a,:) = [ac balAcc(y(te), classify2(abl{a,2}(te,:), abl{a,3}(te,:), P))];
    fprintf('  %-46s %6.1f %6.1f %%\n', abl{a,1}, 100*abA(a,:));
end
[Te, accTrE] = train2(Fm(tr,9:11), rm(tr), y(tr));   yE = classify2(Fm(te,9:11), rm(te), Te);
[Tr, accTrR] = train2(Fm(tr,12:14), rm(tr), y(tr));  yR = classify2(Fm(te,12:14), rm(te), Tr);
[baE, CE] = balAcc(y(te), yE);  [baR, CR] = balAcc(y(te), yR);
fprintf('Trained rule (envelope): outer if R_o > %.1f dB, inner if R_i > %.1f dB (larger margin wins),\n', Te(1:2));
fprintf('  else ball if R_b > %.1f dB or RMS > %.3f g, else healthy\n', Te(3:4));
fprintf('Balanced accuracy  TRAIN: envelope %.1f %%, raw FFT %.1f %%\n', accTrE*100, accTrR*100);
fprintf('Balanced accuracy  TEST : envelope %.1f %%, raw FFT %.1f %%\n', baE*100, baR*100);
fprintf('Per-class recall (test, envelope): '); fprintf('%s %.0f%%  ', cls{1}, 100*CE(1,1)/sum(CE(1,:)));
for c = 2:4, fprintf('%s %.0f%%  ', cls{c}, 100*CE(c,c)/sum(CE(c,:))); end
fprintf('\nFault detection (healthy vs any fault, test): %.1f %% correct\n', 100*mean((y(te)>1) == (yE>1)));
fprintf('Overall accuracy (test): %.1f %% of %d segments\n\n', 100*mean(y(te) == yE), sum(te));
fig = newFig();
subplot(1,2,1); plotCM(CR, cls, sprintf('Raw FFT: %.1f%% balanced accuracy', baR*100));
subplot(1,2,2); plotCM(CE, cls, sprintf('Envelope: %.1f%% balanced accuracy', baE*100));
exportgraphics(fig, fullfile(figdir,'fig3_confusion.png'), 'Resolution', 150);

%% 5) SEVERITY AND NOISE ROBUSTNESS -----------------------------------------
fprintf('==== PHASE 5: SEVERITY AND NOISE ====\n');
fprintf('  Class        size(in)  kurtosis  own indicator (dB)\n');
fig = newFig();  col = lines(4);
for c = 2:4
    sv = zeros(3,2);
    for s = 1:3
        k = y == c & Fm(:,3) == sz(s);  sv(s,:) = [mean(Fm(k,7)) mean(Fm(k,8+[2 3 1]*(c==[2 3 4])'))];
        fprintf('  %-11s  0.%03d  %8.1f  %8.1f\n', cls{c}, sz(s), sv(s,1), sv(s,2));
    end
    subplot(1,2,1); plot(sz/1000, sv(:,1), '-o', 'Color', col(c,:), 'LineWidth', 1.5); hold on;
    subplot(1,2,2); plot(sz/1000, sv(:,2), '-o', 'Color', col(c,:), 'LineWidth', 1.5); hold on;
end
subplot(1,2,1); grid on; xlabel('Fault diameter (in)'); ylabel('Kurtosis'); legend(cls(2:4)); title('Kurtosis vs fault size');
subplot(1,2,2); grid on; xlabel('Fault diameter (in)'); ylabel('Own-class indicator (dB)'); title('Envelope indicator vs fault size');
exportgraphics(fig, fullfile(figdir,'fig4_severity.png'), 'Resolution', 150);
snr = [Inf 20 10 5 0];  accN = zeros(size(snr));
idx = find(te);
for q = 1:numel(snr)
    R = zeros(numel(idx),3);  rN = zeros(numel(idx),1);
    for j = 1:numel(idx)
        seg = Xs(idx(j),:)';
        if isfinite(snr(q)), seg = seg + sqrt(mean(seg.^2)/10^(snr(q)/10))*randn(size(seg)); end
        [fe,Ae] = envSpectrum(seg, fs);  R(j,:) = lpInd(fe, Ae, bearingFreqs(Fm(idx(j),5)));  rN(j) = rmsN(seg);
    end
    accN(q) = balAcc(y(te), classify2(R, rN, Te));
    fprintf('  Added noise SNR %5s dB -> balanced accuracy %.1f %%\n', num2str(snr(q)), accN(q)*100);
end
fig = newFig();  snrPlot = snr;  snrPlot(1) = 30;
plot(snrPlot, accN*100, '-o', 'LineWidth', 1.8);  grid on;  set(gca,'XDir','reverse');
xticks(fliplr(snrPlot));  xticklabels(fliplr({'none','20','10','5','0'}));
xlabel('Added white noise: SNR (dB)'); ylabel('Balanced accuracy (%)'); ylim([0 100]);
title('Robustness of the envelope method to measurement noise');
exportgraphics(fig, fullfile(figdir,'fig5_noise.png'), 'Resolution', 150);

%% 6) SUMMARY ---------------------------------------------------------------
fprintf('\n==== PHASE 6: SUMMARY ====\n');
fprintf('Envelope analysis: %.1f %% balanced accuracy on unseen loads (raw FFT %.1f %%, v1 %.1f %%)\n', baE*100, baR*100, abA(1,2)*100);
fprintf('Healthy-vs-fault detection on unseen loads: %.1f %%\n', 100*mean((y(te)>1) == (yE>1)));
save('RMF_results.mat','Fm','T','Te','Tr','CE','CR','baE','baR','accN','snr','cls','abA','abl');
fprintf('Saved RMF_results.mat and figures/*.png\n');

%% ======================= LOCAL FUNCTIONS ================================
function k = kurtosisN(x), x = x - mean(x); k = mean(x.^4)/mean(x.^2)^2; end
function r = rmsN(x), r = sqrt(mean(x.^2)); end

function R = indicators(f, A, Bf)
% v1 indicator (kept for the ablation): mean amplitude of the first 3 harmonics (+/-2 %) over the 5-1000 Hz median floor, in dB
    fl = median(A(f > 5 & f < 1000));  ff = [Bf.BPFO Bf.BPFI Bf.BSF];  R = zeros(1,3);
    for j = 1:3
        s = 0;
        for h = 1:3
            m = f > h*ff(j)*0.98 & f < h*ff(j)*1.02;  s = s + max(A(m));
        end
        R(j) = 20*log10(s/3/fl);
    end
end

function plotCM(C, cls, ttl)
    P = C./max(sum(C,2),1);  imagesc(P, [0 1]);  colormap(gca, flipud(gray)*0.85 + 0.15);
    axis square;  set(gca,'XTick',1:4,'XTickLabel',cls,'YTick',1:4,'YTickLabel',cls);
    xlabel('Predicted'); ylabel('True'); title(ttl);
    for i = 1:4, for j = 1:4
        text(j,i,sprintf('%d',C(i,j)),'HorizontalAlignment','center','FontWeight','bold', ...
             'Color', [1 1 1]*(P(i,j) < 0.5));
    end, end
end

function gscatter2(x, y, g, cls)
    mk = {'o','s','^','d'};  col = lines(4);  hold on;
    for c = 1:4
        k = g == c;  scatter(x(k), y(k), 28, col(c,:), mk{c}, 'filled', 'MarkerFaceAlpha', 0.6);
    end
    legend(cls, 'Location', 'best');
end

function fig = newFig()
    fig = figure('Color','w','Position',[60 60 1100 700]);
    try, theme(fig,'light'); end %#ok<TRYNC>
end
