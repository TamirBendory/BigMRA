clear; close all; clc;

%% parameters

L = 10; % length of signal
x = sign(randn(L,1)); %x = x - mean(x);
hx = fft(x); hx = hx./abs(hx); x = real(ifft(hx));

N = 10000;
k = 500;

%% generating data

sigma = 0;
y = sigma*randn(N,1);
ind = randi(N, k, 1);
ind  = sort(ind);
indn = ind(1);
for i = 1:length(ind)-1
    if ind(i+1) - ind(i) > 2*L
        indn = [indn; ind(i+1)];
    end
end

ind = indn;

for i = 1:length(ind)
    y( ind(i): ind(i)+L-1) = y( ind(i): ind(i)+L-1) +  x;
end


%%  

Seg = 50;
ys = zeros(Seg,N/Seg);
Ys1 = zeros(Seg,Seg,N/Seg/2);
Ys2 = Ys1;

for i = 1:N/Seg
    
    ys(:,i) = y((i-1)*Seg+1:(i-1)*Seg+Seg);
    
end

for i = 1:N/Seg/2
    Ys1(:,:,i) = ys(:,i)*ys(:,i)';
    Ys2(:,:,i) = ys(:,i+N/Seg/2)*ys(:,i+N/Seg/2)';
end
    
YS1 = mean(Ys1,3);
YS2 = mean(Ys2,3);
YS = YS1/YS2; 

[U V] = eig(YS);    
xest = U(:,1);
figure; hold on; plot(xest); plot(x);

    