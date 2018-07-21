function D = analytic_WignerD(maxL, eu_angles)
eu_angles = sym(eu_angles);
maxL = sym(maxL);

D = cell(maxL+1, 1);

for j = 0:maxL
    D{j+1} = zeros(double(2*j+1), double(2*j+1), 'sym');
    for m = -j:j
        for k = -j:j
            n_min = max(0, k-m);
            n_max = min(j-m, j+k);
            for n = n_min:n_max
                w_jmkn = sqrt(factorial(sym(j+m))*factorial(sym(j-m))*factorial(sym(j+k))*factorial(sym(j-k)));
                w_jmkn = w_jmkn./(factorial(sym(j-m-n))*factorial(sym(j+k-n))*factorial(sym(n+m-k))*factorial(sym(n)));
                D{j+1}(m+j+1,k+j+1) = D{j+1}(m+j+1,k+j+1) + (-1)^n*w_jmkn*cos(eu_angles(2)/2)^(2*j+k-m-2*n)*(-sin(eu_angles(2)/2))^(m-k+2*n);
            end
            D{j+1}(m+j+1, k+j+1) = exp(-1i*m*eu_angles(1))*D{j+1}(m+j+1, k+j+1)*exp(-1i*k*eu_angles(3));
        end
    end
end

% D = cellfun(@(x) double(x), D, 'UniformOutput', 0);