function [M1, M2, M3] = moments_from_harmonics_rand_projs(trials, a_lmi, Psilms_2D, jball_2D, info)
% Evaluate the three first moments of the dataset on a L x L cartesian grid

M1 = 0; M2 = zeros(2*info.L0-1); M3 = zeros(2*info.L0-1,2*info.L0-1,2*info.L0-1,2*info.L0-1);
for t = 1:trials
    proj = gen_rand_proj(a_lmi, Psilms_2D, jball_2D, info);
    [m1, m2, m3] = moments_from_micrograph(proj, info.L0);
    M1 = M1 + m1; M2 = M2 + m2; M3 = M3 + m3;
end
M1 = M1/trials; M2 = M2/trials; M3 = M3/trials;