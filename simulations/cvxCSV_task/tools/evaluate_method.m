function [SNR_out, SDR_out, ISR_out] = evaluate_method(S, I, soi, w)
       
    % SNR_all = zeros(size(w, 2),1);
    % ISR_all = zeros(size(w, 2),1);
    % SDR_all = zeros(size(w, 2),1);

    % for idx = 1:size(w_all,2)
    s_sep = w'*S;
    i_sep = w'*I;
    SNR_out = (mean(s_sep.*conj(s_sep),'all') / mean(i_sep.*conj(i_sep),'all'));
    ISR_out = (mean(i_sep.*conj(i_sep),'all') / mean(s_sep.*conj(s_sep),'all'));
    SDR_out = SDR(s_sep, soi, 0);
    % end
end