function [w, numit] = run_method(X, ini_vals, soi_a, soi_y, method, params)
    splt = strsplit(method,'_');
    method_sel = splt(1);

    switch lower(method_sel)        
        case 'full'
            T = params.T;
            maxIt = params.maxit;
            mu = params.mu;
            a_1_ini = ini_vals.a_1_ini;
            a_T_ini = ini_vals.a_T_ini;
            if ~isfield(params,'update_type')
                [w, ~, ~, numit] = cvxogice_full(X, mu, a_1_ini, a_T_ini, soi_a, soi_y, "maxit", maxIt,'T',T,'nonln',params.nonln);
            else
                [w, ~, ~, numit] = cvxogice_full(X, mu, a_1_ini, a_T_ini, soi_a, soi_y, "maxit", maxIt,'T',T,'nonln',params.nonln, ...
                    'update_type',params.update_type);
            end
        case 'bogice'
            T = params.T;
            maxIt = params.maxit;
            mu = params.mu;
            w_ini = ini_vals.w_ini;
            Nb = size(X,2)/T;

            [w, ~, ~, numit] = bogive_w(X, w_ini, Nb, mu, params.nonln , maxIt, eye(size(X,1)), soi_a, soi_y);
        case 'fastdiva'
            switch params.initype
                case 'a'
                    lambda_t = permute(linspace(0, 1, params.T).',[4, 3, 2, 1]);
                    params.ini = ini_vals.a_1_ini.*(1-lambda_t) + ini_vals.a_T_ini.*(lambda_t);
                case 'w'
                    params.ini = ini_vals.w_ini;
            end
            [w, ~, ~, numit, ~] = fastdiva_pm(X, params);
            w = w(:,:,:,1);
        otherwise
            error('Unknown method')
    end


end