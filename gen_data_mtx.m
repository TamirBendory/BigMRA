function y_mat = gen_data_mtx(y,W,window_shift)

% rearranging the data as a matrix
% input:
% y - (noisy) data
% W - window length
% widnow_shift. If widnow_shift==1, we shift the window at one entry each time 
% output:
% y_mat - the data arranged as a matrix. Each column is a windowed signal

N = length(y); y_stretch = [y ; y(1:W-1)];
Nw = ceil(N/window_shift); % need to verfift the ceiling for non-integer fraction 
y_mat = zeros(W,Nw);

for i = 0:Nw-1
    y_mat(:,i+1) = y_stretch(i*window_shift+1:i*window_shift+W);
end

end
