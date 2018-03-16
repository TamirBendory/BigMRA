%% Setup some fake long signal to compute moments
n = 1e9;
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
%     M3b = M3b + prod(y(k-vals1)); % this is very slow
end
toc

%% Same in parallel -- this is much slower than a regular for-loop on both
%  my laptop and on Latte with 72 workers.
% tic;
M3c = 0;
% parfor k = (1+max(vals1)) : (n+min(vals1))
%     M3c = M3c + y(k)*y(k-shift1)*y(k-shift2);
% % 
% %     suby = y([k k k]-vals1);  % failed attempt to avoid communication overhead
% %     M3c = M3c + suby(1)*suby(2)*suby(3);
% % 
% %     M3c = M3c + y(k)*y(k-7)*y(k+300); % it's actually faster (on both laptop and latte) if shifts are fixed, but still comm. overhead and slower overall.
% end
% toc

%% Parallel again: trying a different way with dummy variables
tic;
M3d = 0;
y1 = y;
y2 = y;
y3 = y;
% On Latte, for n = 1e7, this is faster than the other parfor, but still
% slower than a simple loop. For n = 1e8, it's about the same. For n = 1e9,
% it is slower again. Yeah, so this really doesn't work out...
parfor k = (1+max(vals1)) : (n+min(vals1))
    M3d = M3d + y1(k)*y2(k-shift1)*y3(k-shift2);
end
toc

%% Again, with different dummies, out of despair: this is indeed slower.
% tic;
M3e = 0;
% y1 = y(range1);
% y2 = y(range1-shift1);
% y3 = y(range1-shift2);
% 
% parfor k = 1:numel(range1)
%     M3e = M3e + y1(k)*y2(k)*y3(k);
% end
% toc

%%
fprintf('Relative difference in computed value: %g\n', abs((M3b - M3a) / M3a));
fprintf('Relative difference in computed value: %g\n', abs((M3c - M3a) / M3a));
fprintf('Relative difference in computed value: %g\n', abs((M3d - M3a) / M3a));
fprintf('Relative difference in computed value: %g\n', abs((M3e - M3a) / M3a));