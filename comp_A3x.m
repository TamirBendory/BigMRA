function a = comp_A3x(x,W)
% computing the 3-order AC of x

a = zeros(W);

for i = 0:W-1
    for j = 0:W-1
   
        a(i+1,j+1) = sum(x.*circshift(x,i).*circshift(x,j));
        
    end
end
    