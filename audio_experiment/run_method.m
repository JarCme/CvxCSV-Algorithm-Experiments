function [w] = run_method(X, w_ini, g_1, g_2, method_struct)

switch method_struct.name
    case 'bogive'
        [w, ~, ~, ~, ~] = bogive_w(X, squeeze(w_ini), size(X,2)/method_struct.params.T, method_struct.params.mu, method_struct.params.nonln, method_struct.params.maxit);
    case 'fastdiva'
        params = method_struct.params;
        params.initype = 'w';
        params.ini = w_ini;
        [Wt, ~, ~, ~] = fastdiva(X, params);
        w = Wt(:,:,:,1);
    case 'cvxcsv'  
        [w, ~, ~, ~, ~, ~] = cvxogice_full_ive(X, method_struct.params.mu, ...
                                                    g_1, g_2, 0.5, 0, ...
                                                    "nonln", method_struct.params.nonln, ...
                                                    'T', method_struct.params.T, ...
                                                    'maxit', method_struct.params.maxit);
    case 'overiva'
        [~, W_aux] = overiva_mtlb_wrap(X, "n_iter", method_struct.params.maxit, "W0", w_ini, 'n_src', 1);
        w = permute(W_aux, [2,1,3]);
end

