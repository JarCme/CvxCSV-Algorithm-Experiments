function [w] = randvecangle(v, deg)
    n = norm(v);
    v = v / n;
    d = length(v);
    U = eye(d) - v*v';
%     U = normc(U);
    for c = 1:size(U,2)
        U(:,c) = U(:,c)/norm(U(:,c));
    end
    u = U * randn(d, 1);
    u = u / norm(u);
    w = cos(deg / 180 * pi) * v + sin(deg / 180 * pi) * u;
    w = w * n;
end 

