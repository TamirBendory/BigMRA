clear all;
close all;
clc;
load M3test.mat;
n3 = size(list3, 1);
ns = zeros(n3, 1);

% These two symmetry rules take as input the indices to a 3rd order moment
% entry, and return the indices to a 3rd order moment which has the same
% value (though the indices are usually different.)
% Applying rule 1 and rule 2 in alternance usually produces a cycle of
% length 6, though the cycle can also be shorter (period 1, 2, 3 or 6),
% though I don't understand how the cycles of length 3 can appear. Perhaps
% they are not cycles, but rather lists of length 6 with 3 elements removed
% because of indices out of bounds. AH: no it's fine: the cycle has length
% 6 indeed, but there are only 3 distinct elements in it. Need proper
% language here. First block is invariant under rule 2, and rule 1 points
% to second block. There, rule 2 points to third block. That last one is
% invariant under rule 1.
rule1 = @(a) a([3 4 1 2]);
rule2 = @(a) [-a(1) -a(2) a(3)-a(1) a(4)-a(2)];

g = zeros(n3, 1);

for k = 1 : n3
    
    % Number of repetitions of this value
    n = numel(find(abs(M3 - M3(k)) < 1e-8));
    ns(k) = n;
    
    % Try to find the repetition places
    a = zeros(6, 4);
    a(1, :) = list3(k, :);
    a(2, :) = rule1(a(1, :));
    a(3, :) = rule2(a(2, :));
    a(4, :) = rule1(a(3, :));
    a(5, :) = rule2(a(4, :));
    a(6, :) = rule1(a(5, :));
    % get rid of forbidden lines (there is no wrapping around here)
    keep = zeros(6, 1);
    for line = 1 : 6
        keep(line) = all(a(line, :) >= -top & a(line, :) <= top);
    end
    a = a(logical(keep), :);
    a = unique(a, 'rows');
%     disp(a);
%     disp(n);
    b = list3(abs(M3 - M3(k)) < 1e-8, :);
    b = unique(b, 'rows');
%     disp(b);
    
    if all(size(a) == size(b))
        g(k) = norm(a-b);
    else
        g(k) = Inf;
    end
%     pause;
    
%     if n == 6
%         list3(abs(M3 - M3(k)) < 1e-8, :)
%         break;
%     end
    
end
% hist(ns, 1:15);

fprintf('If this is zero, we had perfect prediction of identical entries: %g\n', norm(g));

