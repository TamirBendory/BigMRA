
addpath('../aspire')
initpath
addpath('~/EMDB_maps')
mol = '1qlq_crop';
K = 2e5

vol = ReadMRC([mol '.mrc']);
parpool('local', maxNumCompThreads)
projections = cryo_project(vol, rand_rots(K));
save(['/scratch/network/eitanl/' mol '_projs_' num2str(K) '.mat'], 'projections', '-v7.3')
