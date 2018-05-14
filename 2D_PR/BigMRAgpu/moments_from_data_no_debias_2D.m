function [M1, M2, M3] = moments_from_data_no_debias_2D(Y, list2, list3)
% Y is an image (a matrix)
% list2 has size n2 x 2
% list3 has size n3 x 4
% Each row of list2 and list3 contains integers (shifts)

if gpuDeviceCount>0 
    	gpu_flag = 1;
else
	gpu_flag = 0;
end

    [N1, N2] = size(Y);

    M1 = sum(Y(:));
    
    n2 = size(list2, 1);

    if gpu_flag
    M2 = zeros(n2, 1,'gpuArray');
    Y = gpuArray(Y);
    else
    M2 = zeros(n2, 1);
    end
    
parfor k = 1 : n2

       %for  k = 1 : n2
        
	shift1 = list2(k, :);
        
        vals1 = [0, shift1(1)];
        range1 = (1+max(vals1)) : (N1+min(vals1));
        vals2 = [0, shift1(2)];
        range2 = (1+max(vals2)) : (N2+min(vals2));
        
        X1 = Y(range1, range2); %#ok<PFBNS>
        X2 = Y(range1-shift1(1), range2-shift1(2));
        
        M2(k) = sum(X1(:) .* X2(:)); % bsxfun
        %T = bsxfun(@times,X1,X2);
        %M2(k) = sum(T(:)); % bsxfun
    end
    
    n3 = size(list3, 1);
    M3 = zeros(n3, 1,'gpuArray');
    
   parfor k = 1 : n3
       %for k = 1 : n3
 
        shifts = list3(k, :);
        shift1 = shifts([1, 2]);
        shift2 = shifts([3, 4]);

        vals1 = [0, shift1(1), shift2(1)];
        range1 = (1+max(vals1)) : (N1+min(vals1));
        vals2 = [0, shift1(2), shift2(2)];
        range2 = (1+max(vals2)) : (N2+min(vals2));
        
        X1 = Y(range1, range2); %#ok<PFBNS>
        X2 = Y(range1-shift1(1), range2-shift1(2));
        X3 = Y(range1-shift2(1), range2-shift2(2));
        
        M3(k) = sum(X1(:) .* X2(:) .* X3(:)); %bsxfun
        %T = bsxfun(@times,X1,bsxfun(@times,X2,X3));
        %M3(k) = sum(T(:)); %bsxfun
    end
      
    if gpu_flag
    M1 = gather(M1);
    M2 = gather(M2);
    M3 = gather(M3);
    end
    
end
