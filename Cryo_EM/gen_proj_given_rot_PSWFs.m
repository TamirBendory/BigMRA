function [I, I_pad, R] = gen_proj_given_rot_PSWFs(Rot, a_lms, psi_lNs, L)

R = RN2RL( getSHrotMtx( Rot, length(a_lms)-1, 'complex' ) );
a_lms_rot = cellfun(@mtimes, a_lms, R, 'UniformOutput', 0);

I = zeros(L);
maxL = length(a_lms)-1;
for l = 0:maxL
    for N = -l:l
        I(:) = I(:) + psi_lNs{l+1}{N+l+1}*a_lms_rot{l+1}(:,N+l+1);
    end
end
I = real(I);

I_pad = zeros(3*L-2);
I_pad(L:2*L-1,L:2*L-1) = I;