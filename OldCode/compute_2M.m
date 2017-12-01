function M = compute_2M(y_mat,ac_flag)

N = size(y_mat,2);

if ac_flag
    
    y1 = y_mat(1,:);
    y1 = repmat(y_mat(1,:),size(y_mat,1),1);
    y = y1.*y_mat;
    M = sum(y,2);
    
%     M = zeros(size(y_mat,1),1);
%     for i = 1:N 
%     M = M + y_mat(1,i)*y_mat(:,i);
%     end
        
else
    
    M = zeros(size(y_mat,1));
    for i = 1:N   %currently with a loop...
        M =  M +  y_mat(:,i)*y_mat(:,i)';
    end
    
end
M = M/N;