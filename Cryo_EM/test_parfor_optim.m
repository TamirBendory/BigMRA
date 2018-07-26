A = 10;
m = 100000; 
n = 10;

a = zeros(m, n, 'double');
tic
for j = 1:n
    for i = 1:m
        a(i,j) = max(abs(eig(rand(A))));
    end
end
toc

a = zeros(m, n, 'double');
tic
parfor j = 1:n
    for i = 1:m
        a(i,j) = max(abs(eig(rand(A))));
    end
end
toc

a = zeros(m, n, 'double');
tic
for j = 1:n
    parfor i = 1:m
        a(i,j) = max(abs(eig(rand(A))));
    end
end
toc

a = zeros(m, n, 'double');
tic
parfor idx = 1:m*n
    a(idx) = max(abs(eig(rand(A))));
end
toc