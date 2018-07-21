function R = get_WignerD_EasySpin(maxL, Rot)

R = cell(maxL+1, 1);
R_angle = eulang(Rot);
for l = 0:maxL
    R{l+1} = wignerd(l, R_angle, '+');
end