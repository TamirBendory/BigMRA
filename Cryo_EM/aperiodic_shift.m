function I_shift = aperiodic_shift(I, shift)
I_shift = zeros(size(I));
sz = size(I);

if sum(shift > sz(1:2)) > 0, return; end

if shift(1) >= 0 && shift(2) >= 0
    I_shift(shift(1)+1:end, shift(2)+1:end, :) = I(1:sz(1)-shift(1), 1:sz(2) - shift(2), :);
elseif shift(1) >= 0 && shift(2) < 0
    I_shift(shift(1)+1:end, 1:sz(2)-abs(shift(2)), :) = I(1:sz(1)-shift(1), abs(shift(2))+1:end, :);
elseif shift(1) < 0 && shift(2) >= 0
    I_shift(1:sz(1)-abs(shift(1)), shift(2)+1:end, :) = I(1+abs(shift(1)):end, 1:sz(2)-shift(2), :);
else % shift(1) < 0 && shift(2) < 0
    I_shift(1:sz(1)-abs(shift(1)), 1:sz(2)-abs(shift(2)), :) = I(1+abs(shift(1)):end, 1+abs(shift(2)):end, :);
end