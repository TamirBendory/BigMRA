function M2 = Second_moment_from_data_stupid_implementation(Y,W,sigma)
% Y is an image (a matrix)

   [N1, N2] = size(Y);
    M2 = zeros(W);

    for shift1 = -(W-1)/2:(W-1)/2
        for shift2 = -(W-1)/2:(W-1)/2

        M2(shift1+(W+1)/2,shift2+(W+1)/2) = sum(sum(Y.*circshift(Y,[shift1,shift2])));

        end
    end

     M2((W+1)/2,(W+1)/2) = M2((W+1)/2,(W+1)/2) -sigma^2*N1*N2;
end

