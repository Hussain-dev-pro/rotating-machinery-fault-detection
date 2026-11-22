function [P, acc] = train2(R, Dt, y)
%TRAIN2  Thresholds for classify2 by coordinate descent (balanced accuracy).
nd = size(Dt,2);  X = [R Dt];
P = [median(R,1) max(Dt,[],1)+1];  acc = -1;
for it = 1:6
    for p = 1:numel(P)
        g = linspace(min(X(:,p)), max(X(:,p)), 120);
        for v = g
            Q = P;  Q(p) = v;
            a = balAcc(y, classify2(R, Dt, Q));
            if a > acc + 1e-12, acc = a; P = Q; end
        end
    end
end
end
