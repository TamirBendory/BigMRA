        function main
        startnow;



        n=10
        x = randn(n,1)


        c = zeros(n,n)


        k=2
        l=4

        kl_del = abs(k - l)


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
