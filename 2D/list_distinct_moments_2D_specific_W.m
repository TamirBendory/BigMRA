% Script to generate list2 and list 3 for specific W

W = 5;
[list2, list3] = list_distinct_moments_2D(W);
 %save(sprintf('lists_W_%d.mat', W), 'list2', 'list3', 'W');