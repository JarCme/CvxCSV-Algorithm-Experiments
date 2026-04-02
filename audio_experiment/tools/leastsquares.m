function alpha = leastsquares(yd,y)
    M = size(yd,3);
    alpha = zeros(M,1);
    for k = 1:M
        alpha(k) = (yd(:,:,k)*y(:,:,k)')/(y(:,:,k)*y(:,:,k)');
    end
end

