function y = p_pdf(s, a, y) % complex generalized gaussian distribution (CGGD) pdf
%     a = 0.5; % laplace - pro jednoduchost zafixovano 
%     y = 0; % circular
    rho = gamma(2/a)/gamma(1/a);
    num = a .* rho .* exp(- ( (rho/2)/(y.^2-1) .* (y.*s.^2 + y.*(conj(s).^2) - 2.*s.*conj(s) )).^a);
    den = pi .* gamma( 1./a ) .* ( 1- y.^2)^(1/2);
    y = num/den;
end