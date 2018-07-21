function [a_aligned, R_rec, cost] = align_vol_coeffs(a_lms, a_lms_ref, trials, par_flag)

problem.M = rotationsfactory(3);
problem.cost = @(Rxyz) cost_rot_alignment(a_lms, a_lms_ref, Rxyz);

R_init = get_init_random_trials(problem, trials, par_flag);
opts.tolgradnorm = 1e-3;
[R_rec, cost] = trustregions(problem, R_init, opts);

J = diag([1;1;-1]);

problem.cost = @(Rxyz) cost_rot_alignment(a_lms, a_lms_ref, J*Rxyz);
R_init = get_init_random_trials(problem, trials, par_flag);
[R_rec_ref, cost_ref] = trustregions(problem, R_init, opts);
if cost_ref < cost
    R_rec = J*R_rec_ref;
    cost = cost_ref;
end

R = RN2RL(getSHrotMtx(R_rec, length(a_lms)-1, 'complex'));
a_aligned = cellfun(@mtimes, a_lms, R, 'UniformOutput', 0);