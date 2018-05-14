function [M1, M2, M3] = moments_from_data_no_debias_2D_v3(Y, list2, list3)
% Y is an image (a matrix)
% list2 has size n2 x 2
% list3 has size n3 x 4
% Each row of list2 and list3 contains integers (shifts)

n2 = size(list2, 1);
if isempty(list3)
    M3 = [];
else
    n3 = size(list3, 1);
end

gpu_flag = (gpuDeviceCount>0);

if gpu_flag > 0
    M2 = zeros(n2, 1,'gpuArray');
    %    M3 = zeros(n3, 1,'gpuArray');
else
    M2 = zeros(n2, 1);
    %   M3 = zeros(n3, 1);
end

[N1, N2] = size(Y);

if gpu_flag
    Y = gpuArray(Y);
    %X1 = zeros(N1,N2,n2,'gpuArray');
    %X2 = zeros(N1,N2,n2,'gpuArray');
end
    %X1 = zeros(N1,N2,n2);
    %X2 = zeros(N1,N2,n2);
%end

M1 = sum(Y(:));

%parfor k = 1 : n2
for k = 1 : n2
    
    shift1 = list2(k, :);
    
    vals1 = [0, shift1(1)];
    range1 = (1+max(vals1)) : (N1+min(vals1));
    vals2 = [0, shift1(2)];
    range2 = (1+max(vals2)) : (N2+min(vals2));
    
    %L1 = length(range1);
    %L2 = length(range2);
    
    X1 = Y(range1, range2); 
    X2 = Y(range1-shift1(1), range2-shift1(2));
    %X1(1:L1,1:L2,k) = Y(range1, range2); %#ok<PFBNS>
    %X2(1:L1,1:L2,k) = Y(range1-shift1(1), range2-shift1(2));
 
    %X1 = X1(:);
    %X2 = X2(:);
    %M2(k) = sum(bsxfun(@times,X1,X2));
    M2(k) = sum(X1(:) .* X2(:)); % bsxfun
    %T = bsxfun(@times,X1,X2);
    %M2(k) = sum(T(:)); % bsxfun
end

%M2 = sum(bsxfun(@times,reshape(X1,N1*N2,[]),reshape(X2,N1*N2,[])),1);
% 
% if 0
%     
%     if gpu_flag
%         X1 = zeros(N1,N2,n3,'gpuArray');
%         X2 = zeros(N1,N2,n3,'gpuArray');
%         X3 = zeros(N1,N2,n3,'gpuArray');
%     else
%         X1 = zeros(N1,N2,n3);
%         X2 = zeros(N1,N2,n3);
%         X3 = zeros(N1,N2,n3);
%     end
%     
    %parfor k = 1 : n3
    for k = 1 : n3
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
        
        %X1 = X1(:);
        %X2 = X2(:);
        %X3 = X3(:);
        
        M3(k) = sum(X1(:) .* X2(:) .* X3(:)); %bsxfun
        %T = bsxfun(@times,X1,bsxfun(@times,X2,X3));
        %M3(k) = sum(bsxfun(@times,X1,bsxfun(@times,X2,X3))); %bsxfun
    end
%end
% 
% M3 = sum(bsxfun(@times,reshape(X3,N1*N2,[]),bsxfun(@times,reshape(X1,N1*N2,[]),reshape(X2,N1*N2,[])),1));
% 
% 
% if gpu_flag > 0
%     M1 = gather(M1);
%     M2 = gather(M2);
%     M3 = gather(M3);
% end

if gpu_flag > 0
    M1 = gather(M1);
    M2 = gather(M2);
    M3 = gather(M3);
end
end
