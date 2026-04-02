function [w, a1_all, aT_all, it, w_all, Shat] = cvxogice_full_ive(X, mu, a1_ini, aT_ini, soi_a, soi_y, options)
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
    end
    %
    [d,N,K] = size(X);
    T = options.T;
    Nb = N/T;
    if( (Nb - fix(Nb))~=0 )
        error('N is is not integer multiple of Nb');
    end
    nonln = char(options.nonln);
    X_in = X;
    X = zeros(d,Nb,K,T);
    for t = 1:T
        n_start = (t-1)*Nb+1;
        n_stop  = (t)*Nb;

        X(:,:,:,t) = X_in(:,n_start:n_stop,:);
    end


    lambda_t = permute(linspace(0, 1, T).',[3, 2, 1]);
    % lambda_t = repmat(lambda_t,1,1,1);
    maxit = options.maxit; 

    Cx_t = pagemtimes(X, 'none', X, 'ctranspose')/Nb;
    Cx = mean(Cx_t, 4);
    invCx_all = pageinv(Cx);

    a1_all = a1_ini;%./norm(a1_ini);
    aT_all = aT_ini;%./norm(a1_ini);

    % save c_vals and w from all iterations
    w_all = zeros(d,1,K,maxit+1);
    w_ini = zeros(d,1,K);
    for k = 1:K
        w_ini(:,:,k)  = lcmvweights([a1_all(:,:,k), aT_all(:,:,k)], [1;1], Cx(:,:,k));
    end
    w = w_ini;
    w_all(:,:,:,1) = w;
   

%     crit = 1000;
    for it = 1:maxit
        % tic;
        % clc;
        % disp(it);
        aux = pagemtimes(Cx_t, w);
   
        sigma2_t = pagemtimes(w,'ctranspose',aux,'none');
    
        sigma_t = sqrt(sigma2_t);
   
        shat = pagemtimes(w, 'ctranspose', X, 'none')./sigma_t;
    
        nu = real(mean( psi(shat,nonln, soi_a, soi_y).* shat , 2));
    
        old_grad_all = mean( - (mean((psi(shat, nonln, soi_a, soi_y).*X), 2)./sigma_t)./nu,4);
        
        lambdaa_t = permute(lambda_t,[4,1,2,3]);
        a_t_all = a1_all.*(1-lambdaa_t) + aT_all.*lambdaa_t;
        g_t_all = a_t_all(2:end,:,:,:);
        gamma_t_all = a_t_all(1,:,:,:);
        B_all = cat(2, g_t_all , -gamma_t_all.*repmat(eye(d-1), 1, 1, K, T) );
        zhat_all = pagemtimes(B_all, X);
        Cz_all = pagemtimes(zhat_all, 'none', zhat_all, 'ctranspose')/Nb;
        zhat_all = permute(zhat_all,[1,5,2,3,4]);
        Cz_all = permute(Cz_all,[1,2,5,3,4]);
        inv_Cz_all = pageinv(Cz_all);
        aux_t3_all = pagemtimes(inv_Cz_all,zhat_all);
        X_vec = permute(X,[1,5,2,3,4]);
        t3_tmp_all = cat(1, ...
                permute(mean(pagemtimes(-X_vec(2:end,:,:,:,:),'ctranspose',aux_t3_all,'none'),3),[1,2,4,5,3]), ...
                permute(mean(pagemtimes(X_vec(1    ,:,:,:,:),'ctranspose',aux_t3_all,'none'),3),[1,2,4,5,3]));
        grad_a1_term3_all = mean(t3_tmp_all.*(1-lambdaa_t),4);
        grad_aT_term3_all = mean(t3_tmp_all.*(lambdaa_t),4);

        grada1_term4_all = (d-2)*cat(1,mean((1./conj(gamma_t_all)).*(1-lambdaa_t),4),zeros(d-1,1,K));
        gradaT_term4_all = (d-2)*cat(1,mean((1./conj(gamma_t_all)).*lambdaa_t,4),zeros(d-1,1,K));
        

        aTH_C_aT_all = pagemtimes(pagemtimes(aT_all,'ctranspose',invCx_all,'none'),aT_all);
        aTH_C_a1_all = pagemtimes(pagemtimes(aT_all,'ctranspose',invCx_all,'none'),a1_all);
        a1H_C_aT_all = pagemtimes(pagemtimes(a1_all,'ctranspose',invCx_all,'none'),aT_all);
        a1H_C_a1_all = pagemtimes(pagemtimes(a1_all,'ctranspose',invCx_all,'none'),a1_all);

        A_all = pagemtimes(invCx_all,a1_all).*aTH_C_aT_all - pagemtimes(invCx_all,a1_all).*a1H_C_aT_all - ...
                pagemtimes(invCx_all,aT_all).*aTH_C_a1_all + pagemtimes(invCx_all,aT_all).*a1H_C_a1_all;

        aux_det_all = a1H_C_a1_all.*aTH_C_aT_all - a1H_C_aT_all.*aTH_C_a1_all;

        graddeta1_a_all = -(aTH_C_aT_all.*pagemtimes(invCx_all,a1_all) - aTH_C_a1_all.*pagemtimes(invCx_all,aT_all))./(aux_det_all).^2;
        graddetaT_a_all = -(a1H_C_a1_all.*pagemtimes(invCx_all,aT_all) - a1H_C_aT_all.*pagemtimes(invCx_all,a1_all))./(aux_det_all).^2;

        invCxa1 = pagemtimes(invCx_all,a1_all);
        invCxaT = pagemtimes(invCx_all,aT_all);
        
        aux_grad_a1_all      = (-   pagemtimes(invCxa1,'none',invCxaT,'transpose') ...
                            +   pagemtimes(invCxaT,'none',invCxa1,'transpose'))./aux_det_all + pagemtimes(A_all,'none',graddeta1_a_all,'transpose');

        aux_grad_conj_a1_all = (    conj(invCx_all).*aTH_C_aT_all ...
                            -   conj(invCx_all).*conj(a1H_C_aT_all) ...
                            -   pagemtimes(conj(invCxaT),'none',invCxaT,'transpose') ...
                            +   pagemtimes(conj(invCxaT),'none',invCxa1,'transpose'))./aux_det_all + pagemtimes(conj(A_all),'none',graddeta1_a_all,'transpose');
        
        aux_grad_aT_all      = (-   pagemtimes(invCxaT,'none',invCxa1,'transpose') ...
                            +   pagemtimes(invCxa1,'none',invCxaT,'transpose'))./aux_det_all + pagemtimes(A_all,'none',graddetaT_a_all,'transpose');

        aux_grad_conj_aT_all = (    conj(invCx_all).*a1H_C_a1_all ...
                            -   conj(invCx_all).*conj(aTH_C_a1_all) ...
                            -   pagemtimes(conj(invCxa1),'none',invCxa1,'transpose') ...
                            +   pagemtimes(conj(invCxa1),'none',invCxaT,'transpose'))./aux_det_all + pagemtimes(conj(A_all),'none',graddetaT_a_all,'transpose');

        grad_a_1_all = permute(pagemtimes(old_grad_all,'ctranspose',aux_grad_a1_all,'none') + pagemtimes(old_grad_all,'transpose',aux_grad_conj_a1_all,'none'),[2,1,3]) - grad_a1_term3_all + grada1_term4_all;
        grad_a_T_all = permute(pagemtimes(old_grad_all,'ctranspose',aux_grad_aT_all,'none') + pagemtimes(old_grad_all,'transpose',aux_grad_conj_aT_all,'none'),[2,1,3]) - grad_aT_term3_all + gradaT_term4_all;
        
        a1_all = a1_all + mu*grad_a_1_all;     
        aT_all = aT_all + mu*grad_a_T_all;

        for k = 1:K
            w(:,:,k) = lcmvweights([a1_all(:,:,k), aT_all(:,:,k)], [1;1], Cx(:,:,k));
            w_all(:,:,k,it+1) = w(:,:,k);
        end

    end
    a_t_all = a1_all.*(1-lambdaa_t) + aT_all.*lambdaa_t;
    Shat_tmp = pagemtimes(w,'ctranspose',X,'none');
    Shat_tmp = a_t_all(1,:,:,:).*Shat_tmp;

    Shat = zeros(1,N,K);
    for t = 1:T
        n_start = (t-1)*Nb+1;
        n_stop  = (t)*Nb;

        Shat(:,n_start:n_stop,:) = Shat_tmp(:,:,:,t);
    end
end

function y_out = psi(s, type, a, y) 

    switch lower(type)
        case 'cggd' % CGGD score function
            y_out = cggd(s, a, y);
        case 'sign'
            sp2 = s.*conj(s);
            aux = 1./sqrt(sum(sp2, 3));
            y_out = conj(s).*aux;
        case 'rati'
            sp2 = s.*conj(s);
            aux = 1./(1+sum(sp2,3));
            y_out = conj(s).*aux;
        case 'sign_ice'
            sp2 = s.*conj(s);
            aux = 1./sqrt(sp2);
            y_out = conj(s).*aux;
    end
end


function y_out = cggd(s, a, y) % CGGD score function
    rho = gamma(2/a)/gamma(1/a);
    fterm = (((2*a)*(rho/2).^a)/((y^2-1)^a));
    sterm = ((y*s.^2+y*(conj(s).^2)-2*(s.*conj(s))).^(a-1)) .* (y*s - conj(s));
    y_out = fterm*sterm;
end

