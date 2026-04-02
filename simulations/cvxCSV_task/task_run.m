function store_struct = task_run(exp_settings, methods_list, methods_params) % main(PBS_ARRAY_INDEX, N, d, T, ini_pert, trials, mu, soi_y, tau, methods)

    addpath('methods');
    addpath('tools');

    store_struct = struct();
    store_struct.IDrange = [exp_settings(1).ID, exp_settings(end).ID];
    store_struct.methods_list = ["ini", methods_list];
    store_struct.res = struct();

    for exp_idx = 1:length(exp_settings) 
        tic

        N           = exp_settings(exp_idx).N; 
        d           = exp_settings(exp_idx).d;
        T           = exp_settings(exp_idx).T;
        m_type      = exp_settings(exp_idx).m_type;
        ini_pert    = exp_settings(exp_idx).ini_pert;
        soi_a       = exp_settings(exp_idx).soi_a;
        soi_y       = exp_settings(exp_idx).soi_y;
        tau         = exp_settings(exp_idx).tau;
        non_soi_a   = exp_settings(exp_idx).non_soi_a;
        non_soi_y   = exp_settings(exp_idx).non_soi_y;
        target_SIR  = exp_settings(exp_idx).target_SIR;

        rng_seed    = exp_settings(exp_idx).seed;
        ID          = exp_settings(exp_idx).ID;
        

        rng(rng_seed);
        
        [X, S, I, soi, ini_vals, ~] = prepare_data(d, T, N, m_type.name, m_type.angle, ini_pert, soi_a, soi_y, tau, non_soi_a, non_soi_y, target_SIR);

        [SNR_ini, SDR_ini, ISR_ini] = evaluate_method(S, I, soi, ini_vals.w_ini);
        
        methods_res = struct();

        methods_res(1).SNR = SNR_ini;
        methods_res(1).ISR = ISR_ini;
        methods_res(1).SDR = SDR_ini;
        methods_res(1).NumIt = 1;

        for m_idx = 1:length(methods_list)
            method = methods_list(m_idx);

            [w, NumIt]                  = run_method(X, ini_vals, soi_a, soi_y, method, methods_params(method));
            [SNR_out, SDR_out, ISR_out] = evaluate_method(S, I, soi, w);
        
            methods_res(m_idx+1).SNR = SNR_out;
            methods_res(m_idx+1).ISR = ISR_out;
            methods_res(m_idx+1).SDR = SDR_out;
            methods_res(m_idx+1).NumIt = NumIt;
        end
        store_struct.res(exp_idx).eval_metrics = methods_res;
        store_struct.res(exp_idx).sens_SIR = (sum(abs(S).^2,[2,3])./sum(abs(I).^2,[2,3]));
        store_struct.res(exp_idx).ID = ID;
        disp(toc)

    end
end







