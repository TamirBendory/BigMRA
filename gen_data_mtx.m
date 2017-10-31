function y_mat = gen_data_mtx(y,window_size)

% rearranging the data as a matrix
% input:
% y - (noisy) dat
% window_size
% output:
% y_mat - the data arranged as a matrix

N = length(y);
y_stretch = [y ; y(1:window_size-1)];
Nw = N/window_size; y_mat = zeros(window_size,Nw);

for i = 0:Nw-1
    y_mat(:,i+1) = y_stretch(i*window_size+1:(i+1)*window_size);
end

end
