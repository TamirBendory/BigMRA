        function main
        startnow;



        len=10
        x = randn(len,1)

%
%        limiting formula for second moment
%
        c = zeros(len,len)

        for k=1:len
%
        for l=1:len
%
        kl_del = abs(k - l)

        for i=1:len - kl_del
%
        c(k,l) = c(k,l) + x(i)*x(i+kl_del);
    end

    end
    end


%
%        check if it's circulant matrix of zero-padded x, squared
%
        x2 = zeros(2*len,1);
        x2(1:len) = x


        chk0 = norm(c - c')


        cx2 = circ_mat(x2,2*len)

        cc = cx2*cx2'
        chk0 = norm(c - cc(1:len,1:len))


        stopnow;
        end
%
%
%
%
%
        function cx = circ_mat(x,m)
%
        cx = zeros(m,m);
        for i=1:m
%
        cx(:,i) = circshift(x,i-1);
    end


        end
%
%
%
%
%
        function startnow
%
        delete out13
        diary('out13')
        diary on
%
        format short E
%%%        format long E


%
        rng('default');

        end
%
%
%
%
%
        function stopnow
%
        diary off
        stop

        end
