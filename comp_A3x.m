function a = comp_A3x(x,W)
% computing the 3-order AC of x
% I am using the definition from our paper 

a = zeros(2*W+1);

for i = -(W-1):W-1
    for j = -(W-1):W-1
   
        a(i+W,j+W) = sum(x.*circshift(x,i).*circshift(x,-j));
        
    end
end
    