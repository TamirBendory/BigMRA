function [a_lms, vol_rec] = expand_vol_3D_PSWF(vol, Psi_lms, pts_in_ball)

vol_ft = cfftn(vol);
L = size(vol, 1);

maxL = length(Psi_lms)-1;
vol_rec = zeros(size(vol_ft));

% sz = cellfun(@size, Psi_lms, 'UniformOutput', 0);
% sz{1}(3) = 1;
% 
% Psi_lms = cellfun(@(x) x(:,:), Psi_lms, 'UniformOutput', 0);
% Psi_lms = horzcat(Psi_lms{:});
% 
% a_lms = Psi_lms\vol_ft(pts_in_ball);
% vol_rec(pts_in_ball) = Psi_lms*a_lms;
% 
% a_lms = mat2cell(a_lms, cellfun(@(x)x(2)*x(3), sz), 1);
% a_lms = cellfun(@(x,y) reshape(x, y(2), y(3)), a_lms, sz, 'UniformOutput', 0);

a_lms = cell(maxL+1, 1);
for l = 0:maxL
    a_lms{l+1} = (Psi_lms{l+1}(:,:)'*Psi_lms{l+1}(:,:))\(Psi_lms{l+1}(:,:)'*vol_ft(pts_in_ball));
    vol_rec(pts_in_ball) = vol_rec(pts_in_ball) + Psi_lms{l+1}(:,:)*a_lms{l+1};
    a_lms{l+1} = reshape(a_lms{l+1}, size(Psi_lms{l+1},2), size(Psi_lms{l+1},3));
end