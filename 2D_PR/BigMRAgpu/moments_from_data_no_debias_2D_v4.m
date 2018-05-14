function [M1, M2, M3] = moments_from_data_no_debias_2D_v4(Y, list2, list3)
% Y is an image (a matrix)
% list2 has size n2 x 2
% list3 has size n3 x 4
% Each row of list2 and list3 contains integers (shifts)

n2 = size(list2, 1);
n3 = size(list3, 1);

gpu_flag = (gpuDeviceCount>0);

if gpu_flag > 0
    M2 = zeros(n2, 1,'gpuArray');
    M3 = zeros(n3, 1,'gpuArray');
else
    M2 = zeros(n2, 1);
    M3 = zeros(n3, 1);
end

[N1, N2] = size(Y);

batch_size = 1;

if gpu_flag
    Y = gpuArray(Y);
    X1 = zeros(N1,N2,batch_size,'gpuArray');
    X2 = zeros(N1,N2,batch_size,'gpuArray');
else
    X1 = zeros(N1,N2,batch_size);
    X2 = zeros(N1,N2,batch_size);
end



%% First moment

M1 = sum(Y(:));

%% Second moment

num_batch = ceil(n2/batch_size);

for i = 1:num_batch
    
    I1 = (i-1)*batch_size+1;
    if n2<i*batch_size
        I2 = n2;
        last_iter = n2 - (i-1)*batch_size;
    else
        I2 = i*batch_size;
        last_iter = batch_size;
    end
    
    for k = 1 : last_iter
        
        shift1 = list2((i-1)*batch_size+k, :);
        vals1 = [0, shift1(1)];
        range1 = (1+max(vals1)) : (N1+min(vals1));
        vals2 = [0, shift1(2)];
        range2 = (1+max(vals2)) : (N2+min(vals2));
        
        L1 = length(range1);
        L2 = length(range2);
        
        X1(1:L1,1:L2,k) = Y(range1, range2);
        X2(1:L1,1:L2,k) = Y(range1-shift1(1), range2-shift1(2));
        
    end
    
    XX1 = reshape(X1,N1*N2,[]);
    XX2 = reshape(X2,N1*N2,[]);
    
    if last_iter<batch_size
        XX1 = XX1(:,1:last_iter);
        XX2 = XX2(:,1:last_iter);
    end
    M2(I1:I2) = sum(XX1.*XX2,1);
    %M2(I1:I2) = sum(bsxfun(@times,XX1,XX2),1);
end


%% Third moment

if gpu_flag
    X1 = zeros(N1,N2,batch_size,'gpuArray');
    X2 = zeros(N1,N2,batch_size,'gpuArray');
    X3 = zeros(N1,N2,batch_size,'gpuArray');
else
    X1 = zeros(N1,N2,batch_size);
    X2 = zeros(N1,N2,batch_size);
    X3 = zeros(N1,N2,batch_size);
end

num_batch = ceil(n3/batch_size);

for i = 1:num_batch
    
    I1 = (i-1)*batch_size+1;
    if n3<i*batch_size
        I2 = n3;
        last_iter = n3 - (i-1)*batch_size;
    else
        I2 = i*batch_size;
        last_iter = batch_size;
    end
    
    for k = 1 : last_iter
        shifts = list3((i-1)*batch_size+k, :);
        shift1 = shifts([1, 2]);
        shift2 = shifts([3, 4]);
        
        vals1 = [0, shift1(1), shift2(1)];
        range1 = (1+max(vals1)) : (N1+min(vals1));
        vals2 = [0, shift1(2), shift2(2)];
        range2 = (1+max(vals2)) : (N2+min(vals2));
        
        L1 = length(range1);
        L2 = length(range2);
        
        X1(1:L1,1:L2,k) = Y(range1, range2);
        X2(1:L1,1:L2,k) = Y(range1-shift1(1), range2-shift1(2));
        X3(1:L1,1:L2,k) = Y(range1-shift2(1), range2-shift2(2));
        
    end
    
    XX1 = reshape(X1,N1*N2,[]);
    XX2 = reshape(X2,N1*N2,[]);
    XX3 = reshape(X3,N1*N2,[]);
    
    if last_iter<batch_size
        XX1 = XX1(:,1:last_iter);
        XX2 = XX2(:,1:last_iter);
        XX3 = XX3(:,1:last_iter);
    end
    
    M3(I1:I2) = sum(XX1.*XX2.*XX3,1);
    %M2(I1:I2) = sum(bsxfun(@times,XX3,bsxfun(@times,XX1,XX2)),1);
end

%% redefining the moments as non-gpu variables

if gpu_flag > 0
    M1 = gather(M1);
    M2 = gather(M2);
    M3 = gather(M3);
end

end
