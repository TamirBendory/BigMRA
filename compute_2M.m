function M = compute_2M(y_mat)

M = zeros(size(y_mat,1));
N = size(y_mat,2);
for i = 1:N   %currently with a loop...
    
 M =  M +  y_mat(:,i)*y_mat(:,i)';
    
end

M = M/N;