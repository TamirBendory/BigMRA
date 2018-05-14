function PS = M2PS_partial(M2,list2,m_total,sigma,N,W)

% This function takes the 2nd moment and arrange it as a power spectrum
% N  - the total number of entries
% sigma - noise level 
% m_total - number of signal's ocuurances 
% M2 - the second moment with only half of the entries

M2_full = zeros(W); 

for i = 1:length(list2)
   M2_full(list2(i,1)+(W+1)/2,list2(i,2)+(W+1)/2) = M2(i,:);
   if all(list2(i,:) == [0,0])
       M2_full(list2(i,1)+(W+1)/2,list2(i,2)+(W+1)/2) = M2_full(list2(i,1)+(W+1)/2,list2(i,2)+(W+1)/2) - sigma^2*N;
   end
end

% symmetry rules
M2_full(1:(W-1)/2,(W+1)/2) = flipud(M2_full((W+1)/2+1:end,(W+1)/2));
M2_full((W+1)/2,1:(W-1)/2) = fliplr(M2_full((W+1)/2,(W+1)/2+1:end));
M2_full(1:(W-1)/2,1:(W-1)/2) = rot90(M2_full((W+1)/2+1:end,(W+1)/2+1:end),2);
M2_full((W+1)/2+1:end,1:(W-1)/2) = rot90(M2_full(1:(W-1)/2,(W+1)/2+1:end),2);
%M2_full = M2_full/m_total;
M2_full = circshift(M2_full,[(W+1)/2,(W+1)/2])/m_total;
%M2((W+1)/2,(W+1)/2) = M2((W+1)/2,(W+1)/2) - sigma^2*N; 
%Mest = circshift(M2,[(W+1)/2,(W+1)/2])/m_total;
%Mest(2:end,2:end) = Mest(2:end,2:end)';
%upper_row = Mest(1,2:end);
%Mest(1,2:end) = Mest(2:end,1);
%Mest(2:end,1) = upper_row;
%Mest(Mest<0) = 0;
%Mest(1,1) = Mest(1,1)- N*sigma^2/m_total;
PS = real(fft2(M2_full));
