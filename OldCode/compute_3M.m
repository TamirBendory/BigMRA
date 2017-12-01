function M = compute_3M(y_mat,ac_flag)

L = size(y_mat,1);
N = size(y_mat,2);

if ac_flag % dummy implementation
    
    M = zeros(L);
    
    for n = 1:N
        
        M = M + y_mat(1,n)*(y_mat(:,n)*y_mat(:,n)');
        
    end
    
else
    
    M = zeros(L,L,L);
    
    for n = 1:N
        
        for i = 1:L
            for j = 1:L
                for k = 1:L
                    M(i,j,k) =  M(i,j,k) +  y_mat(i,n)*y_mat(j,n)*y_mat(k,n);
                end
            end
        end
    end
    
end

M = M/N;