function list = symmetry_rules_3rd_order_2D(c3, top)
% Inputs
%   c3: row vector of size 4x1 containing the indices to a 3rd order moment
%   top: maximum value allowed for indices
%
% Output
%   list: a matrix of size nx4, with n ranging from 1 to 6. Each row
%         contains indices to a 3rd order moment which is identical to c3.

    % These two symmetry rules take as input the indices to a 3rd order
    % moment entry, and return the indices to a 3rd order moment which has
    % the same value (though the indices are usually different.)
    %
    % Notice that rule1 and rule2 are their own repsective inverses.
    rule1 = @(a) a([3 4 1 2]);
    rule2 = @(a) [-a(1) -a(2) a(3)-a(1) a(4)-a(2)];
    
    % Applying rule 1 and rule 2 in alternance usually produces a cycle of
    % length 6, though the pattern can contain fewer than 6 distinct
    % elements (possible numbers are 1, 2, 3 and 6). These (maximum) 6 sets
    % of indices form an exhaustive list of all 3rd order moments that are
    % equivalent to the one indexed by c3.
    list = zeros(6, 4);
    list(1, :) = c3;
    list(2, :) = rule1(list(1, :));
    list(3, :) = rule2(list(2, :));
    list(4, :) = rule1(list(3, :));
    list(5, :) = rule2(list(4, :));
    list(6, :) = rule1(list(5, :));
    
    % We filter out sets of indices that exceed the maximum allowed index,
    % which is -top and top.
    keep = zeros(6, 1);
    for line = 1 : 6
        keep(line) = all(list(line, :) >= -top & list(line, :) <= top);
    end
    list = list(logical(keep), :);
    
    % Avoid repetitions in the list.
    list = unique(list, 'rows');

end
