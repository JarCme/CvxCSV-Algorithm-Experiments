clear;clc;close all;

if (exist('res','dir') ~= 7)
    disp('unzipping precomputed results')
    unzip("res.zip")
    disp('done');
end

methods = ["ini","overiva","fastdiva","cvxcsv","bogive"];
blks_lens = [100,200,400,800,1200,1600,2000];

res_SIR = zeros(length(blks_lens),length(methods));
res_SDR = zeros(length(blks_lens),length(methods));

res_path = "res"+filesep;

for method_idx = 1:length(methods)
    method = methods(method_idx);
    for blks_idx = 1:length(blks_lens)
        blks_len = blks_lens(blks_idx);
        i = 1;
        all_SIR = [];
        all_SDR = [];
        for f = 1:40 % 40 pairs
            fol = res_path+string(f)+filesep+string(blks_len)+filesep+method;
            for blk = 1:numel(dir(fol+filesep+"*.mat"))
                l_ =load(fol+filesep+string(blk)+".mat");
                all_SIR(i) = l_.SIR;
                all_SDR(i) = l_.SDR;
                i = i+1;
            end
        end
        res_SIR(blks_idx,method_idx) = median(10*log10(all_SIR));
        res_SDR(blks_idx,method_idx) = median(10*log10(all_SDR));
    end
end

SNR_picked = permute(res_SIR,[3,2,1]);
SDR_picked = permute(res_SDR,[3,2,1]);

n_table(SNR_picked, SDR_picked, blks_lens, methods);