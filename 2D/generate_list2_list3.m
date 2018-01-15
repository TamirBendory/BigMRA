% Call list_distinct_moments_2D for various values of W
for W = 1 : 2 : 51
    
    fprintf('Working on W = %3d ...', W);
    [list2, list3] = list_distinct_moments_2D(W);
    save(sprintf('lists_W_%d.mat', W), 'list2', 'list3', 'W');
    fprintf(' done.\n');
    
end
