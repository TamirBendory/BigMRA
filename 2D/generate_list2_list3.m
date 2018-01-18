function generate_list2_list3(Ws)
% Call list_distinct_moments_2D for various values of W

    assert(all((Ws+1)/2 == round((Ws+1)/2)), 'Only tested for odd W''s.');

    parfor k = 1 : numel(Ws)
        
        W = Ws(k);

        fprintf('Working on W = %3d ...', W);
        [list2, list3] = list_distinct_moments_2D(W);
        save_lists(W, list2, list3);
        
        fprintf(' done.\n');

    end

end

function save_lists(W, list2, list3) %#ok<INUSD>

    save(sprintf('lists_W_%d.mat', W), 'list2', 'list3', 'W');

end
