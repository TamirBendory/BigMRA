%% Setup some fake long signal to compute moments
n = 1e8;
y = randn(n, 1);

shifts = [7, -300];
shift1 = shifts(1);
shift2 = shifts(2);

vals1 = [0, shift1, shift2];


%% The clean, loop-free code is slow...
tic;
range1 = (1+max(vals1)) : (n+min(vals1));
y1 = y(range1);
y2 = y(range1-shift1);
y3 = y(range1-shift2);
M3a = sum(y1 .* y2 .* y3);
% Combining y1,y2,y3,M3 in one instruction does not speed things up.
toc


%% This loop is actually faster! Gains are higher on my laptop (R2017b) than on Latte (R2015a)
tic;
M3b = 0;
for k = (1+max(vals1)) : (n+min(vals1))
    M3b = M3b + y(k)*y(k-shift1)*y(k-shift2);
end
toc

%%
fprintf('Relative difference in computed value: %g\n', abs((M3b - M3a) / M3a));