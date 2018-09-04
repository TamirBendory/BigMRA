function identifiability_of_gamma_and_sigma_Sep4_2018()
% See notes Sep 4, 2018, NB36
L = 10;
labels = zeros(L^3, 3);
counter = 0;
for k1 = 0 : L-1
    for k2 = k1 : L-1
        for k3 = k2 : L-1
            counter = counter + 1;
            labels(counter, :) = [k1, k2, k3]; % only sorted triplets, to account for symmetry of monomials
        end
    end
end
labels = labels(1:counter, :);
nlabels = size(labels, 1);

indices = @(k) labels(k, :);

% labels(linear_index(k1, k2, k3), :) == sort([k1, k2, k3])
% that is: indices composed with linear_index is sort.
function j = linear_index(k1, k2, k3)
    K = sort([k1, k2, k3]);
    % terrible code but not important
    for j = 1 : nlabels
        if all(labels(j, :) == K)
            return;
        end
    end
    error('should not get here...');
end
for k = 1 : size(labels, 1)
    assert(all(indices(linear_index(k1, k2, k3)) == sort([k1, k2, k3])));
end

% Each third-order autocorrelation entry gives access to a certain linear
% combination of monomials of degree 3. There are L^2 entries.
A = zeros(nlabels, L^2);
counter = 0;
for k1 = 0 : L-1
    for k2 = 0 : L-1
        counter = counter + 1;
        % Setup linear combination corresponding to a_x^3[k1, k2]
        for j = max([0, -k1, -k2]) : (L-1+min([0, -k1, -k2]))
            A(linear_index(j, j+k1, j+k2), counter) = 1;
        end
    end
end

% To identify gamma and sigma, this is what we want to express as a linear
% combination of 3rd order autocorrelation entries:
b = zeros(nlabels, 1);
for k1 = 0 : L-1
    for k2 = 0 : L-1
        b(linear_index(k1, k1, k2)) = L;
    end
end

x = A\b;
fprintf('Make sure this is close to zero: %g\n', norm(A*x-b)/norm(x));

for k1 = 0 : L-1
    for k2 = 0 : L-1
        k = k1*L + k2 + 1;
        if abs(x(k)) >= 1e-12
            fprintf('Take (%d, %d) times %g\n', k1, k2, x(k));
        end
    end
end



end