function [pos, rcv_clk_error] = nonlinearLS(pseudorange, INFO)
    N = size(pseudorange,1);

    solution = zeros(N, 4);
    for k = 1:N
        fun = @(x)estimatePosition(x, pseudorange(k,:)', INFO);
        if k==1
            solution(k,:) = lsqnonlin(fun, zeros(1,4));
        else
            solution(k,:) = lsqnonlin(fun, solution(k-1,:));
        end
    end

    pos = solution(:,1:3);
    rcv_clk_error = solution(:,4);
end

function F = estimatePosition(x, pseudorange, INFO)
    rSV = INFO.rSvECEF;
    Ip = physconst("LightSpeed") * INFO.dtIono;
    T  = physconst("LightSpeed") * INFO.dtTropo; 
    pseudorange_corrected = pseudorange - Ip - T + physconst("LightSpeed")*INFO.dtS;

    rUser = x(1:3);
    rcv_clk_delay = x(4);
    F = vecnorm(rSV - rUser, 2, 2) + rcv_clk_delay - pseudorange_corrected;
end