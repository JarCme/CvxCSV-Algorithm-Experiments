function  out = g_pdf(s, a, y)    
    rho = gamma(2/a)/gamma(1/a);
    fterm = (((2*a)*(rho/2).^a)/((y^2-1)^a));
    sterm = ((y*s.^2+y*(conj(s).^2)-2*(s.*conj(s))).^(a-1)) .* (y*s - conj(s));
    out = fterm*sterm;
end

