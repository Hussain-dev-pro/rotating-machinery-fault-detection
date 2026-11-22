%% RMF_Animation.m  -  animated diagnosis dashboard on unseen test data
% Run RMF_Main.m first (creates RMF_results.mat and data/).
% For one file of each class at 2 hp (a load NOT used for training) it shows:
%   - the vibration signal scrolling past,
%   - the envelope spectrum with the bearing defect frequencies marked,
%   - the three envelope indicators and the overall RMS against their trained
%     thresholds (two-stage rule of RMF_Main) and the verdict.
% Output: animation/RMF_dashboard.gif and a key frame per class.
clear; clc; close all;
S = load('RMF_results.mat');  Te = S.Te;  cls = S.cls;
fs = 12000;  segN = 12000;  files = [99 107 120 132];  rpmNom = 1750;   % 2 hp, 0.007 in
outdir = 'animation';  if ~isfolder(outdir), mkdir(outdir); end
gif = fullfile(outdir,'RMF_dashboard.gif');
fig = figure('Color','w','Position',[40 40 1250 680]);  try, theme(fig,'light'); end %#ok<TRYNC>
axT = axes(fig,'Position',[0.05 0.70 0.62 0.22]);
axE = axes(fig,'Position',[0.05 0.10 0.62 0.48]);
axB = axes(fig,'Position',[0.74 0.36 0.22 0.56]);
hV  = annotation(fig,'textbox',[0.70 0.06 0.29 0.22],'EdgeColor','none','FontSize',15,'FontWeight','bold', ...
                 'HorizontalAlignment','center','VerticalAlignment','middle');
nF = 30;  n = 0;
for c = 1:4
    [x, rpm] = loadCWRU(files(c), 'data', fs, rpmNom, c == 1);
    seg = x(1:segN);  Bf = bearingFreqs(rpm/60);
    [fe,Ae,band] = envSpectrum(seg, fs);  R = lpInd(fe, Ae, Bf);  rmsV = sqrt(mean((seg-mean(seg)).^2));
    yh = classify2(R, rmsV, Te);  ok = yh == c;
    % static parts of the frame
    cla(axE); plot(axE, fe, Ae, 'b', 'LineWidth', 1); hold(axE,'on'); grid(axE,'on');
    xlim(axE,[0 500]); ylim(axE,[0 1.15*max(Ae(fe>5 & fe<500))]);
    mk = {Bf.BPFO,[0.85 0.1 0.1],'BPFO'; Bf.BPFI,[0 0.6 0],'BPFI'; Bf.BSF,[0.85 0.45 0],'BSF'};
    for q = 1:3, for h = 1:3
        xline(axE, h*mk{q,1}, ':', 'Color', mk{q,2}, 'LineWidth', 1.4);
    end, end
    text(axE, Bf.BPFO, 0.97*axE.YLim(2), ' BPFO', 'Color', mk{1,2}, 'FontWeight','bold');
    text(axE, Bf.BPFI, 0.90*axE.YLim(2), ' BPFI', 'Color', mk{2,2}, 'FontWeight','bold');
    text(axE, Bf.BSF,  0.83*axE.YLim(2), ' BSF',  'Color', mk{3,2}, 'FontWeight','bold');
    xlabel(axE,'Envelope frequency (Hz)'); ylabel(axE,'Envelope amplitude (g)');
    title(axE, sprintf('Envelope spectrum  (band %.1f-%.1f kHz chosen by kurtosis)', band/1e3));
    cla(axB); bar(axB, 1:3, R, 0.6, 'FaceColor', [0.3 0.45 0.7]); hold(axB,'on');
    for q = 1:3, plot(axB, q+[-0.4 0.4], Te(q)*[1 1], 'r-', 'LineWidth', 2.5); end
    set(axB,'XTick',1:3,'XTickLabel',{'Outer','Inner','Ball'}); ylabel(axB,'Indicator (dB above local floor)');
    title(axB,{'Indicators vs thresholds (red)', sprintf('RMS %.3f g  (alarm %.3f g)', rmsV, Te(4))}); grid(axB,'on'); ylim(axB,[0 max([R Te])*1.25]);
    vc = [0.1 0.55 0.15]*ok + [0.8 0.1 0.1]*(~ok);
    set(hV, 'String', {sprintf('TRUE: %s', cls{c}), sprintf('DIAGNOSIS: %s', cls{yh}), ...
        sprintf('file %d.mat, 2 hp, %.0f rpm', files(c), rpm)}, 'Color', vc);
    for k = 1:nF
        i0 = round((k-1)/nF*(segN-1200)) + 1;  t = (i0:i0+1199)/fs*1e3;
        plot(axT, t, seg(i0:i0+1199), 'k');  grid(axT,'on');  xlim(axT, t([1 end]));
        ylim(axT, 1.1*max(abs(seg))*[-1 1]);  xlabel(axT,'Time (ms)');  ylabel(axT,'Accel. (g)');
        title(axT, sprintf('Drive-end vibration, unseen 2 hp test data  (%s)', cls{c}));
        drawnow;  n = n + 1;  fr = getframe(fig);  img = fr.cdata;
        [A,map] = rgb2ind(img(1:2:end,1:2:end,:), 128, 'nodither');
        if n == 1, imwrite(A,map,gif,'gif','LoopCount',Inf,'DelayTime',0.12);
        else,      imwrite(A,map,gif,'gif','WriteMode','append','DelayTime',0.12); end
        if k == nF, imwrite(img, fullfile(outdir, sprintf('frame_%d_%s.png', c, strrep(cls{c},' ','_')))); end
    end
end
fprintf('Saved %s (%d frames)\n', gif, n);

