        function main
        startnow;



        len=10
        x = randn(len,1)


        c = zeros(len,len)


        k=2
        l=4

        kl_del = abs(k - l)

        for i=1:len - kl_del
%
        c(k,l) = c(k,l) + x(i)*x(i+kl_del);
    end

        c(k,l)

        stopnow;
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
