function [f, g, store] = costgrad_micro_moments_projs...
    (a_lms, L_list, L, Rots, Psilms_2D, jball_2D,PSWF_quad_int, points_inside_the_circle, m3_micro, gamma, store)

a_lms = accumarray(L_list+1, a_lms, [max(L_list)+1, 1], @(x) {x});
for l = 0:length(a_lms)-1
    a_lms{l+1} = reshape(a_lms{l+1}, [], 2*l+1);
end

if ~isfield(store, 'dm')
    [~, ~, m3_harms] = moments_from_harmonics_projs(Rots, a_lms, Psilms_2D, jball_2D, L, PSWF_quad_int, points_inside_the_circle);
    dm = cellfun(@(x,y) gamma*x(1:size(y,1), 1:size(y,2)) - y, m3_harms(1:length(m3_micro)), m3_micro, 'UniformOutput', 0);
    dm = cellfun(@(x) x(:), dm, 'UniformOutput', 0);
    store.dm = dm;
else
    dm = store.dm;
end


if ~isfield(store, 'f')
    f = sum(cellfun(@(x) norm(x(:))^2, dm));
    store.f = f;
else
    f = store.f;
end

if nargout >= 2
    if ~isfield(store, 'g')
        tmp = cellfun(@(x,y) x(:,1:length(y))*y, G(1:length(dm)), dm, 'UniformOutput', 0);
        g = 0;
        for k = 1:length(tmp)
            g = g + tmp{k};
        end
        store.g = g;
    else
        g = store.g;
    end
end