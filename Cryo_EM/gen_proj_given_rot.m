function [I, R] = gen_proj_given_rot(Rot, a_lms, Psilms_2D, jball_2D, L)

R = RN2RL( getSHrotMtx( Rot, length(a_lms)-1, 'complex' ) );
a_lms_rot = cellfun(@mtimes, a_lms, R, 'UniformOutput', 0);
a_lms_rot = cellfun(@(x) x(:,1:2:end), a_lms_rot, 'UniformOutput', 0);

I = recover_from_ALM_v4_given_Psilms_2D(a_lms_rot, floor(L/2), jball_2D, Psilms_2D, length(a_lms)-1, L);
I = real(icfft2(I));