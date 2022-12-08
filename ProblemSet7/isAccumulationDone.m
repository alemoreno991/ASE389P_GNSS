function flag = isAccumulationDone(tauj, Ta)
    persistent k
    if isempty(k)
        k = 1;
    end

    if tauj > k*Ta
        k = k + 1;
        flag = true;
    else
        flag = false;
    end
end
