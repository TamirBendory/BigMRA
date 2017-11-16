clear; close; clc;

load('errXP1');

k_length = 10;
k_vec = round(logspace(2,5,k_length));
e = mean(err,2);
ln = 1.5;

figure; plot(k_vec,e,'linewidth',ln);
%set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
axis tight
xlabel('K')
ylabel('error')

saveas(gcf,'errXP1')
pdf_print_code(gcf, 'errXP1', 11)