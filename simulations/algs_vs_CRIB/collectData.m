function [ISR_picked, SNR_picked, SDR_picked, params_to_get_vals, methods] = collectData(folder_name, param_to_get)
    path_ = [char(folder_name),filesep];
    path = [path_,'exp_settings.mat'];
    load(path, 'exp_settings');
    load(path, 'methods');
    load(path, 'methods_params');
    
    dir_res = dir(path_);
    files_name = string({dir_res.name});
    files_name = files_name(3:end-2);
    batch_length = length(exp_settings)/length(files_name);
    res = [];
    for f_idx = 1:length(files_name)
        batch_file_name = ['ID_',num2str((f_idx-1)*batch_length+1),'-',num2str((f_idx)*batch_length),'.mat'];
        load([path_,batch_file_name],'store_struct')
        res = [res, store_struct.res];
    end
    
    ISR = [];
    SNR = [];
    SDR = [];
    get_param_array = [];
    for exp_idx = 1:length(exp_settings)
        ISR = [ISR;res(exp_idx).eval_metrics.ISR];
        SNR = [SNR;res(exp_idx).eval_metrics.SNR];
        SDR = [SDR;res(exp_idx).eval_metrics.SDR];
        if exp_settings(exp_idx).ID == exp_idx
            get_param_array = [get_param_array; getfield(exp_settings(exp_idx),param_to_get)];
        else 
            error('Something went wrong');
        end
    end
    params_to_get_vals = sort(unique(get_param_array));
    ISR_picked = zeros(size(ISR,1)/length(params_to_get_vals), size(ISR,2), length(params_to_get_vals));
    SNR_picked = zeros(size(SNR,1)/length(params_to_get_vals), size(SNR,2), length(params_to_get_vals));
    SDR_picked = zeros(size(SDR,1)/length(params_to_get_vals), size(SDR,2), length(params_to_get_vals));
    for get_vals_idx = 1:length(params_to_get_vals)
        ISR_picked(:,:,get_vals_idx) = ISR(get_param_array == params_to_get_vals(get_vals_idx),:);
        SNR_picked(:,:,get_vals_idx) = SNR(get_param_array == params_to_get_vals(get_vals_idx),:);
        SDR_picked(:,:,get_vals_idx) = SDR(get_param_array == params_to_get_vals(get_vals_idx),:);
    end
    methods = ["ini_lcmp",methods];
end