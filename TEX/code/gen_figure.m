% Code for generating the XC figure

rng(1953);

L = 11;
N = 100;
ind = [15, 55];
x = ones(L,1);
y = zeros(N,1);
y(ind(1):ind(1)+L-1) = x;
y(ind(2):ind(2)+L-1) = x;
sigma1 = 0.5;
y1 = y + sigma1*randn(N,1);
sigma2 = 3;
y2 = y + sigma2*randn(N,1);

xc = xcorr(x,y); 
xc1 = xcorr(x,y1); 
xc2 = xcorr(x,y2); 

ln = 1.5;
fs = 14;

figure(321); 

%subplot(321); plot(1:N,y,'linewidth',ln); 
subplot(321); plot(1:80,y(1:80),'linewidth',ln); 
t1 = title('Measurement');
titlePos = get( t1 , 'position');
titlePos(2) = 1.1;
set( t1 , 'position' , titlePos);
ylabel('\sigma = 0'); set(get(gca,'ylabel'),'rotation',0)
set(get(gca,'ylabel'),'Units', 'Normalized','Position',[-.2,0.4,0]);
set(get(gca,'ylabel'),'FontSize',fs);
axis tight

%subplot(322); plot(1:N,xc(L+1:N+L),'linewidth',ln); 
subplot(322); plot(1:80,xc(32:111),'linewidth',ln); 
t2 = title('Cross-correlation with the true signal');
titlePos = get( t2 , 'position');
titlePos(2) = 12;
set( t2 , 'position' , titlePos);
axis tight

%ylim([ylimXlow,ylimXhigh]);

%subplot(323); plot(1:N,y1,'linewidth',ln); 
subplot(323); plot(1:80,y1(1:80),'linewidth',ln); 
%ylim([-10,15]);
ylabel(['\sigma = ' , num2str(sigma1)]); set(get(gca,'ylabel'),'rotation',0)
set(get(gca,'ylabel'),'Units', 'Normalized','Position',[-.2,0.4,0]);
set(get(gca,'ylabel'),'FontSize',fs);
axis tight

%subplot(324); plot(1:N,xc1(L+1:N+L),'linewidth',ln); 
subplot(324); plot(1:80,xc1(32:111),'linewidth',ln); 
axis tight

%ylim([ylimXlow,ylimXhigh]);

%subplot(325); plot(1:N,y2,'linewidth',ln); 
subplot(325); plot(1:80,y2(1:80),'linewidth',ln); 
%ylim([-10,15]);
ylabel(['\sigma = ' , num2str(sigma2)]); set(get(gca,'ylabel'),'rotation',0)
set(get(gca,'ylabel'),'Units', 'Normalized','Position',[-.2,0.4,0]);
set(get(gca,'ylabel'),'FontSize',fs);
axis tight

%subplot(326); plot(1:N,xc2(L+1:N+L),'linewidth',ln); 
subplot(326); plot(1:80,xc2(32:111),'linewidth',ln); 
axis tight
%ylim([ylimXlow,ylimXhigh]);

%pdf_print_code(gcf,'XC_example.pdf', 14)

