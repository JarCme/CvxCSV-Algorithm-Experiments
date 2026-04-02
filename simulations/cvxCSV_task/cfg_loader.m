function [exp_setting, method_list, methods_params] = cfg_loader(cfg_path)

    N = []; d = []; T = []; a1aT_type = []; ini_pert = []; soi_a = []; soi_y = []; tau = []; trials = []; method_list = []; methods_params = []; target_SIR = [];
    non_soi_a = []; non_soi_y = [];

    run(cfg_path)
    var_to_test =   [   isempty(N), isempty(d), isempty(T), ...
                        isempty(ini_pert), isempty(soi_a), ...
                        isempty(tau), isempty(trials), ...
                        isempty(method_list), isempty(methods_params), ...
                        isempty(a1aT_type)];
    if any(var_to_test)
        error('Incomplete config setting');
    end

    if isempty(target_SIR)
        target_SIR = nan;
    end

    if isempty(non_soi_a)
        non_soi_a = 1;
    end
    if isempty(non_soi_y)
        non_soi_y = 0;
    end

       
    [N, d, T, soi_a ,soi_y ,tau] = ndgrid(N, d, T, soi_a, soi_y, tau);
    
    exps_count = length(N)*trials;
    exps_total_idxs = 1:(exps_count*length(ini_pert));

    N = repmat( N(:), trials, 1);
    d = repmat( d(:), trials, 1);
    T = repmat( T(:), trials, 1);
    soi_a = repmat( soi_a(:), trials, 1);
    soi_y = repmat( soi_y(:), trials, 1);
    tau = repmat( tau(:), trials, 1);
    exp_setting = [];
    for ini_idx = 1:length(ini_pert)
        % for m = 1:length(method_list)
            exp_setting_tmp = struct();

            for i = 1:exps_count
                exp_setting_tmp(i,1).N            = N(i);
                exp_setting_tmp(i,1).d            = d(i);
                exp_setting_tmp(i,1).T            = T(i);
                exp_setting_tmp(i,1).soi_a        = soi_a(i);
                exp_setting_tmp(i,1).soi_y        = soi_y(i);
                exp_setting_tmp(i,1).tau          = tau(i);
                exp_setting_tmp(i,1).non_soi_a    = non_soi_a;
                exp_setting_tmp(i,1).non_soi_y    = non_soi_y;
                exp_setting_tmp(i,1).seed         = exps_total_idxs( ((ini_idx-1)*exps_count) + i); % i
                exp_setting_tmp(i,1).m_type       = a1aT_type;
                exp_setting_tmp(i,1).ID           = exps_total_idxs( ((ini_idx-1)*exps_count) + i);
                exp_setting_tmp(i,1).target_SIR   = target_SIR;

                % exp_setting_tmp(i,1).method       = method_list(m);
                exp_setting_tmp(i,1).ini_pert     = ini_pert(ini_idx);
            end
    
            exp_setting = [exp_setting; exp_setting_tmp];
        % end
    end
    
end

