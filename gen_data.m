function [y, yc, ind] = gen_data(x,N,k,L,sigma) 

yc = zeros(N,1);
ind = randi([L,N-L], k, 1);

ind  = sort(ind);
indn = ind(1);
for i = 1:length(ind)-1
    if ind(i+1) - ind(i) > 4*L
        indn = [indn; ind(i+1)];
    end
end

ind = indn;

for i = 1:length(ind)
    yc( ind(i): ind(i)+L-1) = yc( ind(i): ind(i)+L-1) +  x;
end

y = yc + sigma*randn(N,1);

end