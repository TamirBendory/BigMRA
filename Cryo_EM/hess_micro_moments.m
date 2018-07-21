function [Hu, store] = hess_micro_moments(x, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, m1_micro, lambda, u, store)

if ~isfield(store, 'J')
    [~,~,store] = costgrad_micro_moments(x, x_lists, L, B, W, B_lists, M2_quants, m3_micro, m2_micro, m1_micro, lambda, store);
end
J = store.J;

Hu = 2*J.'*(J*u);