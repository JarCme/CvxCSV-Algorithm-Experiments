function [w, a1, aT, it] = cvxogice_full(X, mu, a1_ini, aT_ini, soi_a, soi_y, options)
    arguments
        X double
        mu double
        a1_ini double
        aT_ini double
        soi_a double % pdf "shape", laplace -> soi_a = 0.5
        soi_y double % circularity, circular -> soi_y = 0
        options.nonln = 'cggd';
        options.maxit (1,1) {mustBeNumeric, mustBeInteger, mustBePositive} = 1000;
        options.T (1,1) {mustBeNumeric, mustBeInteger, mustBePositive} = 4;
        options.update_type = 'vanilla';
    end
    %
    [d, N] = size(X);
    T = options.T;
    Nb = N/T;
    if( (Nb - fix(Nb))~=0 )
        error('N is is not integer multiple of Nb');
    end
    nonln = char(options.nonln);
    X = reshape(X, d, Nb, T);

    update_type = options.update_type;
    
    lambda_t = permute(linspace(0, 1, T).',[3, 2, 1]);

    maxit = options.maxit; 

    Cx_t = pagemtimes(X, 'none', X, 'ctranspose')/Nb;
    Cx = mean(Cx_t, 3);
    invCx = inv(Cx);
    X_aux = permute(X,[1,4,2,3]);

    a1 = a1_ini;
    aT = aT_ini;


    w  = lcmvweights([a1, aT], [1;1], Cx);
  
   
    G1 = zeros(size(a1));
    GT = zeros(size(aT));
    for it = 1:maxit
        aux = pagemtimes(Cx_t, w);
   
        sigma2_t = pagemtimes(w,'ctranspose',aux,'none');
    
        sigma_t = sqrt(sigma2_t);

        a_t = a1.*(1-lambda_t) + aT.*lambda_t;
        g_t = a_t(2:end,:,:);
        gamma_t = a_t(1,:,:);

        aTH_C_aT = aT'*invCx*aT;
        aTH_C_a1 = aT'*invCx*a1;
        a1H_C_aT = a1'*invCx*aT;
        a1H_C_a1 = a1'*invCx*a1;
    
        A = invCx*a1*aTH_C_aT- invCx*a1*a1H_C_aT - invCx*aT*aTH_C_a1 + invCx*aT*a1H_C_a1;
    
        aux_det = a1H_C_a1*aTH_C_aT - a1H_C_aT*aTH_C_a1; 
    
        graddeta1_a = -(aTH_C_aT*invCx*a1 - aTH_C_a1*invCx*aT)/(aux_det)^2;
        graddetaT_a = -(a1H_C_a1*invCx*aT - a1H_C_aT*invCx*a1)/(aux_det)^2;
    
        aux_grad_a1      = (-   invCx*a1*(invCx*aT).' ...
                            +   invCx*aT*(invCx*a1).')/aux_det + A*graddeta1_a.';
        aux_grad_conj_a1 = (    conj(invCx)*aTH_C_aT ...
                            -   conj(invCx)*conj(a1H_C_aT) ...
                            -   conj(invCx)*conj(aT)*(aT.'*conj(invCx)) ...
                            +   conj(invCx)*conj(aT)*(invCx*a1).')/aux_det + conj(A)*graddeta1_a.';
    
        aux_grad_aT      = (-   invCx*aT*(invCx*a1).' ...
                            +   invCx*a1*(invCx*aT).')/aux_det + A*graddetaT_a.';
        aux_grad_conj_aT = (    conj(invCx)*a1H_C_a1 ...
                            -   conj(invCx)*conj(aTH_C_a1) ...
                            -   conj(invCx)*conj(a1)*(a1.'*conj(invCx)) ...
                            +   conj(invCx)*conj(a1)*(invCx*aT).')/aux_det + conj(A)*graddetaT_a.';
    
        shat = pagemtimes(w, 'ctranspose', X, 'none')./sigma_t;
    
        nu = real(mean( psi(shat,nonln, soi_a, soi_y).* shat , 2));
    
        B = cat(2, g_t , -gamma_t.*repmat(eye(d-1), 1, 1, T) );

        zhat = pagemtimes(B, X);

        Cz = pagemtimes(zhat, 'none', zhat, 'ctranspose')/Nb;
        Cz = permute(Cz, [1,2,4,3]);
        zhat = permute(zhat, [1,4,2,3]);
        inv_Cz = pageinv(Cz);
        aux_t3 = pagemtimes(inv_Cz,zhat);
        
        grad_a1_term3 = mean(cat(1, ...
                            permute(mean(pagemtimes(-X_aux(2:end,:,:,:),'ctranspose',aux_t3,'none'),3),[1,2,4,3]),...
                            permute(mean(pagemtimes(X_aux(1    ,:,:,:),'ctranspose',aux_t3,'none'),3),[1,2,4,3])).*(1-lambda_t),3);
        grad_aT_term3 = mean(cat(1, ...
                            permute(mean(pagemtimes(-X_aux(2:end,:,:,:),'ctranspose',aux_t3,'none'),3),[1,2,4,3]),...
                            permute(mean(pagemtimes(X_aux(1    ,:,:,:),'ctranspose',aux_t3,'none'),3),[1,2,4,3])).*lambda_t,3);

        grada1_term4 = (d-2)*[mean((1./conj(gamma_t)).*(1-lambda_t),3);zeros(d-1,1)];
        gradaT_term4 = (d-2)*[mean((1./conj(gamma_t)).*lambda_t,3);zeros(d-1,1)];

        old_grad = mean( - (mean((psi(shat, nonln, soi_a, soi_y).*X), 2)./sigma_t)./nu,3);


        grad_a_1 = (old_grad'*aux_grad_a1 + old_grad.'*aux_grad_conj_a1).' - grad_a1_term3 + grada1_term4;
        grad_a_T = (old_grad'*aux_grad_aT + old_grad.'*aux_grad_conj_aT).' - grad_aT_term3 + gradaT_term4;


%-------------------------- Vanila gradient update
        switch char(lower(update_type))
            case 'vanilla'
                % if mod(it,2) == 0
                    a1 = a1 + mu*grad_a_1;
                % else
                    aT = aT + mu*grad_a_T;
                % end
                
            case 'adagrad'
                G1 = G1 + grad_a_1.^2;
                GT = GT + grad_a_T.^2;
    
                a1 = a1 + mu*grad_a_1./sqrt(G1);
                aT = aT + mu*grad_a_T./sqrt(G1);
            otherwise
                error('unknown update type.');
        end      

        w = lcmvweights([a1, aT], [1;1], Cx);

        % crit = norm([grad_a_1;grad_a_T]);

        % if(crit < 1e-6)
        %     break;
        % end
    end
end

% function y = p(s, a, y) % complex generalized gaussian distribution (CGGD) pdf
% %     a = 0.5; % laplace
% %     y = 0; % circular
%     rho = gamma(2/a)/gamma(1/a);
%     num = a .* rho .* exp(- ( (rho/2)/(y.^2-1) .* (y.*s.^2 + y.*(conj(s).^2) - 2.*s.*conj(s) )).^a);
%     den = pi .* gamma( 1./a ) .* ( 1- y.^2)^(1/2);
%     y = num/den;
% end

function y_out = psi(s, type, a, y) % CGGD score function

    switch lower(type)
        case 'cggd'
            y_out = cggd(s, a, y);
        case 'sign'
            sp2 = s.*conj(s);
            aux = 1./sqrt(sum(sp2, 3));
            y_out = conj(s).*aux;
    end
end


function y_out = cggd(s, a, y) % CGGD score function
    rho = gamma(2/a)/gamma(1/a);
    fterm = (((2*a)*(rho/2).^a)/((y^2-1)^a));
    sterm = ((y*s.^2+y*(conj(s).^2)-2*(s.*conj(s))).^(a-1)) .* (y*s - conj(s));
    y_out = fterm*sterm;
end

