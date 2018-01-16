L = 11;
W = L;
shift1 = [-10, 2]; % shift

X = randn(L, L);

% inner product over matrices
inner = @(X, Y) X(:)'*Y(:);
% Zero-padding operator
ZP = @(X) [X zeros(L, W); zeros(W, L+W)];
% Adjoint of the zero-padding operator
ZPadj = @(Y) Y(1:L, 1:L);
% 2D circular shift operator
CS = @(X, k) circshift(X, k);
% Adjoint of the 2D circular shift operator
CSadj = @(X, k) circshift(X, -k);

% Original code
inner(ZP(X), CS(ZP(X), shift1))

% This does a lot of work but worth it if looking at all shifts anyway, and
% it does work for negative k's, in -(L-1) .. (L-1).
F = conv2(X, rot90(X, 2));
F(L+shift1(1), L+shift1(2))

% Explicit formula, works for -(L-1) <= k <= L-1
if shift1(1) >= 0
    range11 = (1+shift1(1)):L;
    range21 = 1:(L-shift1(1));
else
    range11 = 1:(L+shift1(1));
    range21 = (1-shift1(1)):L;
end
if shift1(2) >= 0
    range12 = (1+shift1(2)):L;
    range22 = 1:(L-shift1(2));
else
    range12 = 1:(L+shift1(2));
    range22 = (1-shift1(2)):L;
end
inner(X(range11, range12), X(range21, range22))

% This explicit formula should generalize better
range1 = max(1, 1+shift1(1)) : min(L, L+shift1(1));
range2 = max(1, 1+shift1(2)) : min(L, L+shift1(2));
inner(X(range1, range2), X(range1-shift1(1), range2-shift1(2)))


%% Going for third order
shift1 = [0,  2];
shift2 = [5, -3];
ZPX = ZP(X);
CX1 = CS(ZPX, shift1);
CX2 = CS(ZPX, shift2);
T1 = CX1.*CX2;
inner(ZPadj(T1), X)

[L1, L2] = size(X);
vals1 = [0, shift1(1), shift2(1)];
range1 = (1+max(vals1)) : (L1+min(vals1));
vals2 = [0, shift1(2), shift2(2)];
range2 = (1+max(vals2)) : (L2+min(vals2));
X1 = X(range1, range2);
X2 = X(range1-shift1(1), range2-shift1(2));
X3 = X(range1-shift2(1), range2-shift2(2));
sum(X1(:) .* X2(:) .* X3(:))