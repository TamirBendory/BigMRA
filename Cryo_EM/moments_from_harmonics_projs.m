function [M1, M2, M3] = moments_from_harmonics_projs(Rots, a_lms, Psilms_2D, jball_2D, L, PSWF_Nn_p, PSWF_quad_int)
% Evaluate the three first moments of the dataset on a L x L cartesian grid

trials = size(Rots, 3);
projs = zeros(L, L,  trials); 
for t = 1:trials
    projs(:,:,t) = gen_proj_given_rot(Rots(:,:,t), a_lms, Psilms_2D, jball_2D, L);
end
[M1, M2, M3] = moments_from_micrograph_steerable(projs, L, PSWF_Nn_p, PSWF_quad_int);