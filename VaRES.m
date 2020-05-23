function [VaR,ES] = VaRES(Sample,VaRLevel)

    N = length(Sample);
    % k = ceil(X) rounds each element of X to the nearest integer greater than or equal to that element.
    k = ceil(N*VaRLevel);
    % Sorts the elements of Sample in descending order.
    z = sort(Sample,'descend');
    % VaR est le quantile 1-VaRLevel de la distribution de Sample.
    VaR = z(k);
    
    if k < N
        % Average of the returns beyond the VaR level 
       ES = ((k - N*VaRLevel)*z(k) + sum(z(k+1:N)))/(N*(1 - VaRLevel));
    else
       ES = z(k);
    end

end