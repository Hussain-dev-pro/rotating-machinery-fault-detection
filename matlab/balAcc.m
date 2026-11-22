function [ba, C] = balAcc(y, yh)
%BALACC  Balanced accuracy (mean per-class recall) and 4x4 confusion matrix.
C = accumarray([y(:) yh(:)], 1, [4 4]);
ba = mean(diag(C)./max(sum(C,2),1));
end
