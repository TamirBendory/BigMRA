clear all; close all;
addpath('~/Documents/kam_cryo')
addpath('../SHT')
addpath('../aspire')
initpath
addpath('~/EMDB_maps')
addpath(genpath('../pswf_t_mod'))
% addpath(genpath('../manopt'))

info.mol = '1qlq_crop';
vol = double(ReadMRC([info.mol '.mrc'])); %read volume
L = length(vol);
M = 7420; % size of a micrograph
num_micros = 1000; % number of micrographs
projs_per_micrograph = ceil(0.15*(M/L)^2); % upper bound on # of projection per micrograph

I = zeros(M, M, num_micros, 'double'); % micrographs
gamma = zeros(num_micros, 1, 'double'); % occupancy factors

if isempty(gcp('nocreate')), parpool('local', maxNumCompThreads); end

parfor t = 1:num_micros % for each micrograph
    % Generate projections for current micrograph:
    projs_curr = cryo_project(vol, rand_rots(projs_per_micrograph));

    num_projs = 0; % counter for projections
    inds_possible = cart2inds([M,M], 1:M-L+1, 1:M-L+1); % possible indices for upper-left corner
    curr_micro = zeros(M,M,'double');
    while ~isempty(inds_possible)
        % Put projection at random possible index, increment counter:
        [off_1, off_2] = ind2sub([M,M],inds_possible(randi(length(inds_possible))));
        curr_micro(off_1:off_1+L-1, off_2:off_2+L-1) = projs_curr(:,:,num_projs+1);        
        num_projs = num_projs + 1;
        
        % Delete indices near current projection from possible indices:
        inds_to_del = cart2inds([M,M], ...
            off_1-2*(L-1):off_1+2*(L-1), off_2-2*(L-1):off_2+2*(L-1));
        inds_possible = setdiff(inds_possible, inds_to_del);
    end
    I(:,:,t) = curr_micro;
    gamma(t) = num_projs*L^2/M^2;
end

save(['/scratch/network/eitanl/' info.mol '_micros_' num2str(num_micros) '_M_' num2str(M) '.mat'], '-v7.3')
