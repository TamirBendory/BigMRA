function f = cost_rot_alignment(a_lms, a_lms_ref, Rxyz)

R_l = RN2RL(getSHrotMtx(Rxyz, length(a_lms)-1, 'complex'));
a_lms_rot = cellfun(@mtimes, a_lms, R_l, 'UniformOutput', 0);

f = sum(cellfun(@(x,y) norm(x(:) - y(:))^2, a_lms_rot, a_lms_ref));