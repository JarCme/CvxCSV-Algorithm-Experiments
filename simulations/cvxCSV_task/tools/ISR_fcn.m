function [ISR_dB_BICE, ISR_dB_CSV, ISR_dB_CSV_nonstat, ISR_dB_linearCSV] = ISR_fcn(N,T,d,Cz,soi_a,soi_y,tau)
    Nb = N/T;
    sigma_soi = (tau*ones(T,1) + (1-tau)*sin((1:T).'./T*pi./2));
    K = ((soi_a^2)*gamma(2/soi_a))/((1-(soi_y^2))*(gamma(1/soi_a)^2));
    
    %% unit variance of soi over all blocks
    ISR_dB_BICE = 10*log10((T*(d-1))/(N*(K-1)) );
    % ISR_dB_CMV = 10*log10(  ((d-1)/N) * (1/(K-1) + (T-1) / K));
    ISR_dB_CSV = 10*log10( (d-1)/(N*(K-1)) );
    
    %% varying variance of soi across blocks
    T_CSV_tmp = zeros(d-1); 
    for t = 1:T
        T_CSV_tmp = T_CSV_tmp + (1/sigma_soi(t).^2)*Cz(:,:,t);
    end
    T_CSV = trace(     1/sum(sigma_soi.^2) * inv(T_CSV_tmp) * sum(Cz,3));
    ISR_dB_CSV_nonstat = 10*log10( T/(N*(K-1)) * T_CSV);
    
    %% varying variance of soi across blocks (linear model)
    
    CRLB_C = -(T/2) * [eye(d-1), eye(d-1)];
    CRLB_B = -(T/2) * [eye(d-1); eye(d-1)];
    CRLB_D = zeros(d-1);
    CRLB_A = zeros(2*(d-1),2*(d-1));
    
    K_s = (((soi_a^2)*gamma(2/soi_a))/((1-(soi_y^2))*(gamma(1/soi_a)^2)))./(sigma_soi.^2);
    
    for t = 1:T
        CRLB_D = CRLB_D + K_s(t)*Cz(:,:,t);
    
        lambda = (t-1)/(T-1);
        tmp_S = sigma_soi(t).^2 * inv(Cz(:,:,t));
        CRLB_A = CRLB_A + [   (1-lambda).^2 * tmp_S     , (1-lambda)*lambda * tmp_S ;...
                    (1-lambda)*lambda * tmp_S , (lambda).^2 * tmp_S ]; 
    end
    
    CRLB = (1/Nb) * inv( CRLB_D - CRLB_C*inv(CRLB_A)*CRLB_B );
    TMP = zeros(d-1);
    for t = 1:T
        TMP = TMP + Cz(:,:,t)*CRLB;
    end
    
    ISR_dB_linearCSV = 10*log10((1/sum(sigma_soi.^2)) * trace(TMP));
end

