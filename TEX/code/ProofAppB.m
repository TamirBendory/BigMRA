% script to check the proof of Appendix B

L = 40; 
x = randn(L,1);
x = [x ; zeros(L,1)];

% Verifying F1

A00 = sum(x.^3);
Aii = zeros(L-1,1);
A0i = zeros(L-1,1);
 
 for i = 1:L-1
    Aii(i) =  sum(x.*(circshift(x,i).^2));
 end
 for i = 1:L-1
    A0i(i) =  sum((x.^2).*(circshift(x,i)));
 end


% check F1 indentity (should be zero)
sum(x)*norm(x)^2 - (A00+sum(Aii)+sum(A0i))

% Verifying F2

% computing additional terms
Aij = zeros(L-1,L-1);

for j = 1:L-1
    for i = j+1: L-1 
    Aij(i,j) =  sum(x.*circshift(x,i).*circshift(x,j));
    end
end

% check F2 indentity (should be zero)
sum(x)^3 - (A00+3*sum(Aii)+3*sum(A0i)+6*sum(Aij(:)))

