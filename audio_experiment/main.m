clear;clc;
addpath("tools");

addpath("methods"+filesep+"overiva");
addpath("methods"+filesep+"fastdiva");
addpath("methods"+filesep+"cvxcsv");
addpath("methods"+filesep+"bogive");

tic;
fs_target = 16000;
rng(1)
folder_name = 'res';
mkdir(folder_name)
copyfile([mfilename,'.m'], folder_name);

%% parameters

pick_shift = 100;
numframes_start = 10;
numframes_end = 10;
init_pert = 0;

T = 10;
maxIt = 600;


methods(1).name = 'fastdiva';
methods(1).params.approach = 'u'; 
methods(1).params.numsig = 1;
methods(1).params.initype = 'w';
methods(1).params.T = T;
methods(1).params.nonln = 'rati';
methods(1).params.scaling = 'n';
methods(1).params.maxit = maxIt;

methods(end+1).name = 'bogive';
methods(end).params.nonln = 'rati';
methods(end).params.maxit = maxIt;
methods(end).params.T = T;
methods(end).params.mu = 0.1;

methods(end+1).name = 'cvxcsv';
methods(end).params.nonln = 'rati';
methods(end).params.maxit = maxIt;
methods(end).params.T = T;
methods(end).params.mu = 0.001;

% may need to be computed separately
% methods(end+1).name = 'overiva';
% methods(end).params.maxit = maxIt;

l_ = load('pairs_aux_flac.mat'); % result of file_name_gen script
soi_int_pairs = l_.soi_int_pairs; % pairs where interfererer is static and soi is moving


%% data load
for pair_idx = 1:size(soi_int_pairs,1)
    if ~exist(string(folder_name)+filesep+string(num2str(pair_idx)),"dir")
        mkdir(string(folder_name)+filesep+string(num2str(pair_idx)))
    end
    disp(string(pair_idx)+"   "+soi_int_pairs(pair_idx,2)+"   "+soi_int_pairs(pair_idx,1));
    soi_f_name = soi_int_pairs(pair_idx,2);
    int_f_name = soi_int_pairs(pair_idx,1);
    [x, fs] = audioread(soi_f_name);
    
    soi = x((2.1*fs):(end-fs*2.6),:); % crop silence
    soi = resample(soi, fs_target, fs);
   
    [x, fs] = audioread(int_f_name);
    interference = x((2.1*fs):(end-fs*2.6),:);  % crop silence
    interference = 0.7*resample(interference, fs_target, fs);
    
    fs = fs_target;
    x = soi+interference;
    
    nfft = 256;
    
    overlap_len = 3*nfft/4;
    
    win = [zeros(nfft/2,1);hamming(nfft,'periodic')];
    win2 = [hamming(nfft,'periodic');zeros(nfft/2,1)];
    const = max(win+win2);
    
    X_all = stft(x,"FFTLength",nfft,"Window",hamming(nfft,"periodic"),"OverlapLength",overlap_len,"FrequencyRange","onesided");
    S_all = stft(soi,"FFTLength",nfft,"Window",hamming(nfft,"periodic"),"OverlapLength",overlap_len,"FrequencyRange","onesided");
    I_all = stft(interference,"FFTLength",nfft,"Window",hamming(nfft,"periodic"),"OverlapLength",overlap_len,"FrequencyRange","onesided");
    
    for block_len = [100,200,400,800,1200,1600,2000]
        if ~exist(string(folder_name)+filesep+string(num2str(pair_idx))+filesep+string(num2str(block_len)),"dir")
            mkdir(string(folder_name)+filesep+string(num2str(pair_idx))+filesep+string(num2str(block_len)))
        end
        f_start_idx_all = 1:pick_shift:(length(X_all)-block_len);
        f_stop_idx_all = (1:pick_shift:(length(X_all)-block_len))+block_len-1;
        
        SIR_cvx = zeros(1,length(f_stop_idx_all));
        SDR_cvx = zeros(1,length(f_stop_idx_all));
        Shat_cvx_out_blks = zeros(1,pick_shift,size(X_all,1),length(f_start_idx_all));
        Ssep_cvx_out_blks = zeros(1,pick_shift,size(X_all,1),length(f_start_idx_all));
        Isep_cvx_out_blks = zeros(1,pick_shift,size(X_all,1),length(f_start_idx_all));
        
        parfor blk_idx = 1:length(f_start_idx_all)
                % profile on;
            % tic;
            % mkdir(string(folder_name)+filesep+string(num2str(pair_idx))+filesep+string(block_len)+filesep+string(blk_idx))
            f_start_idx = f_start_idx_all(blk_idx);
            f_stop_idx = f_stop_idx_all(blk_idx);
            
            X = X_all(:, f_start_idx:f_stop_idx,:);
            S = S_all(:, f_start_idx:f_stop_idx,:);
            I = I_all(:, f_start_idx:f_stop_idx,:);
            
            X = permute(X, [3,2,1]);
            S = permute(S, [3,2,1]);
            I = permute(I, [3,2,1]);
            
            s_time = istft(permute(S,[3,2,1]),fs,"Window",hamming(nfft,'periodic'),"OverlapLength",overlap_len,"FrequencyRange","onesided","Method","ola",'ConjugateSymmetric',true);
            i_time = istft(permute(I,[3,2,1]),fs,"Window",hamming(nfft,'periodic'),"OverlapLength",overlap_len,"FrequencyRange","onesided","Method","ola",'ConjugateSymmetric',true);
            
            g_1 = ones(size(X,1),1,size(X,3)); 
            g_2 = ones(size(X,1),1,size(X,3));
    
            sS = S;
            PL_1 = mean(squeeze(conj(sS(1,1:numframes_start,:)).*sS(1,1:numframes_start,:)),1);
            PL_2 = mean(squeeze(conj(sS(1,(end-numframes_end+1):end,:)).*sS(1,(end-numframes_end+1):end,:)),1);
            for m = 2:size(X,1)
                PLR_1=mean(squeeze(conj(sS(1,1:numframes_start,:)).*sS(m,1:numframes_start,:)),1);
                PLR_2=mean(squeeze(conj(sS(1,(end-numframes_end+1):end,:)).*sS(m,(end-numframes_end+1):end,:)),1);
                g_1(m,1,:) = (PLR_1./PL_1);
                g_2(m,1,:) = (PLR_2./PL_2);
            end
    
            eps1 = crandn([size(X,1), 1, size(X,3)]);
            eps2 = crandn([size(X,1), 1, size(X,3)]);
            for k = 1:size(X, 3)
                eps1(:,:,k) = eps1(:,:,k)./norm(eps1(:,:,k));
                eps2(:,:,k) = eps2(:,:,k)./norm(eps2(:,:,k));
            end
    
            g_1 = g_1 + eps1*init_pert;
            g_2 = g_2 + eps2*init_pert;
        
            Cx = pagemtimes(X, 'none', X, 'ctranspose')/size(X, 2);
    
            w_ini = zeros(size(X,1),1,size(X,3));
            for k = 1:size(X,3)
                w_ini(:,:,k) = lcmvweights([g_1(:,:,k), g_2(:,:,k)], [1;1], Cx(:,:,k));
            end

            [SIR_ini, SDR_ini] = eval_w(X, S, I, s_time, w_ini, fs, nfft, overlap_len);
            if ~exist(string(folder_name)+filesep+string(num2str(pair_idx))+filesep+string(block_len)+filesep+"ini","dir")
                mkdir(string(folder_name)+filesep+string(num2str(pair_idx))+filesep+string(block_len)+filesep+"ini");
            end
            s1 = struct("SIR",SIR_ini,"SDR",SDR_ini); 
            target_ini_f_name = string(folder_name)+filesep+string(num2str(pair_idx))+filesep+string(block_len)+filesep+"ini"+filesep+string(blk_idx)+".mat";
            if ~exist(target_ini_f_name,'file')
                save(target_ini_f_name,"-fromstruct",s1);
            end
            for m = 1:length(methods)
                target_folder_name = string(folder_name)+filesep+string(num2str(pair_idx))+filesep+string(block_len)+filesep+string(methods(m).name)
                if ~exist(target_folder_name,"dir")
                    mkdir(target_folder_name);
                end
                target_method_f_name = string(folder_name)+filesep+string(num2str(pair_idx))+filesep+string(block_len)+filesep+string(methods(m).name)+filesep+string(blk_idx)+".mat";
                if ~exist(target_method_f_name,'file')
                    w = run_method(X, w_ini, g_1, g_2, methods(m));
                    [SIR, SDR] = eval_w(X, S, I, s_time, w, fs, nfft, overlap_len);
                    
                    s2 = struct("SIR",SIR,"SDR",SDR); 
                    save(target_method_f_name,"-fromstruct",s2);
                end
            end
            
        end
    end
end
