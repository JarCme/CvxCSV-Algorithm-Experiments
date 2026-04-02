function mtlb_fun(PBS_ARRAY_INDEX, PBS_ARRAY_TOTAL, PBS_NCPUS, config_name)
	maxNumCompThreads(PBS_NCPUS);

    config_name = char(config_name);

    config_file_path = ['configs', filesep, config_name];
    output_path = ['outputFiles', filesep, replace(config_name,'.','_(dot)_'), filesep];
    
    if ~exist(output_path,'dir')
        mkdir(output_path);
    end

    [exp_settings, methods, methods_params] = cfg_loader(config_file_path);

    if PBS_ARRAY_INDEX == 1
        copyfile(config_file_path, [output_path,'config.m']);
        save([output_path, 'exp_settings.mat'], "exp_settings", "methods", "methods_params");
    end

    total_experiments = length(exp_settings);
    per_array = total_experiments/PBS_ARRAY_TOTAL;

    exp_indexes = ((PBS_ARRAY_INDEX-1)*per_array+1):((PBS_ARRAY_INDEX)*per_array);

    store_struct = task_run(exp_settings(exp_indexes), methods, methods_params);

    save([output_path, 'ID_', num2str(exp_indexes(1)), '-', num2str(exp_indexes(end)) ,'.mat'], "store_struct");
    
end
