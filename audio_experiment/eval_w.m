function [SIR_out, SDR_out] = eval_w(X, S, I, s_time, w, fs, nfft, overlap_len)

    Shat_ini = pagemtimes(w(:,:,:),'ctranspose',X,'none');
    Ssep_ini = pagemtimes(w(:,:,:),'ctranspose',S,'none');
    Isep_ini = pagemtimes(w(:,:,:),'ctranspose',I,'none');

    alpha = leastsquares(X(1,:,:),Shat_ini); % projection on 1. mic
   
    s_sep   = istft(permute(Ssep_ini.*permute(alpha,[3,2,1]),[3,2,1]),  fs,"Window",hamming(nfft,'periodic'),"OverlapLength",overlap_len,"FrequencyRange","onesided","Method","ola",'ConjugateSymmetric',true);
    i_sep   = istft(permute(Isep_ini.*permute(alpha,[3,2,1]),[3,2,1]),  fs,"Window",hamming(nfft,'periodic'),"OverlapLength",overlap_len,"FrequencyRange","onesided","Method","ola",'ConjugateSymmetric',true);
    
    SIR_out = var(s_sep)./var(i_sep);
    SDR_out = SDR(s_sep,s_time(:,1),0);
end

