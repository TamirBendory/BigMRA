function coeff_cell = PSWF_coeff_convert(coeff, ang_freq)

ang_freq_unique = unique(ang_freq);
coeff_cell = cell(length(ang_freq_unique), 1);

for ii = 1:length(ang_freq_unique)
    coeff_cell{ang_freq_unique(ii) + 1} = coeff(ang_freq == ang_freq_unique(ii),:);
end
