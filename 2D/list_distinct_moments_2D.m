function [list2, list3] = list_distinct_moments_2D(W)

    % Code limited to odd W for now.
    assert((W-1)/2 == round((W-1)/2), 'W assumed odd in this code.');

    % Maximum index (minimum index is -top)
    top = (W-1)/2;
    
    % Build list2 for 2nd order moments: we can generate this one directly.
    list2 = zeros((W^2+1)/2, 2);
    k = 0;
    for k1 = 0 : top
        k = k + 1;
        list2(k, :) = [k1, 0];
    end
    for k1 = (-top) : top
        for k2 = 1 : top
            k = k + 1;
            list2(k, :) = [k1, k2];
        end
    end
    
    % Build list3 for 3rd order moments; for this one, we need to go over
    % all possible moments and check symmetry rules to exclude most.
    list3 = zeros(W^4, 4);
    k = 0;
    for k1 = (-top) : top
        for k2 = (-top) : top
            for l1 = (-top) : top
                for l2 = (-top) : top
                    candidate = [k1, k2, l1, l2];
                    equivalents = symmetry_rules_3rd_order_2D(candidate, top);
                    if isempty(intersect(list3(1:k, :), equivalents, 'rows'))
                        k = k + 1;
                        list3(k, :) = candidate;
                    end
                end
            end
        end
    end
    list3 = list3(1:k, :);

end
