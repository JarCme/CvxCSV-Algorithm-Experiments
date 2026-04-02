function [X, S, I, soi, ini_vals, true_vals] = prepare_data(d, T, N, mixing_type, a1_aT_angle, ini_pert, soi_a, soi_y, tau, non_soi_a, non_soi_y, target_SIR)
    if (T == -1)
        T = N;
    end
    
    Nb = N/T;
    sigma_soi = tau*ones(T,1) + (1-tau)*sin((1:T).'./T*pi./2);
    
    % non_soi_y = 0; % background sources gamma (0 -> circular)
    % non_soi_a = 1; % background sources alpha (1 -> gaussian)

    [A, w_true] = generate_linear_mixing_parameters(d, T, 'complex', mixing_type, a1_aT_angle);

    soi_stationary_non_blks = cggd_rand(soi_a, 1, N, soi_y);
    interferences_non_blks = zeros(d-1, N);
    
    for d_ = 1:d-1
        interferences_non_blks(d_,:) = cggd_rand(non_soi_a, 1, N, non_soi_y);
    end

    % (d,N) -> (d,Nb,T)
    soi_blocks = reshape(soi_stationary_non_blks, 1, Nb, T);
    interferences_blocks = reshape(interferences_non_blks, d-1, Nb, T);

    % nonstationary soi sigma on blocks (according to tau and line 3)
    soi = soi_blocks.*reshape(sigma_soi, 1, 1, T);
    interferences = interferences_blocks;
        
    sources = cat(1, soi, interferences);
    
    if ~isnan(target_SIR)
        sc = sqrt(exp((target_SIR/10)*log(10))/ mean(  sum(abs(pagemtimes(A(:,1,:), soi)).^2,[2,3])./sum(abs(pagemtimes(A(:,2:end,:), interferences)).^2,[2,3])   ));
        A(:,1,:) = sc * A(:,1,:);
    end

    X = pagemtimes(A, sources);
    S = pagemtimes(A(:,1,:), soi);
    I = pagemtimes(A(:,2:end,:), interferences);
    % disp(mean(abs(X-(S+I)).^2,'all')); % -> 0
       
    eps_1 = crandn(d, 1);
    eps_1 = eps_1./norm(eps_1);
    eps_T = crandn(d, 1);
    eps_T = eps_T./norm(eps_T);
    
    a_1_ini = A(:,1,1) + eps_1*(ini_pert + tau);
    a_T_ini = A(:,1,T) + eps_T*(ini_pert + tau);

    Cx = mean(pagemtimes(X, 'none', X, 'ctranspose')/Nb, 3); 
    Cx_true =  pagemtimes(A(:,1,:),'none',A(:,1,:),'ctranspose').*reshape(sigma_soi,1,1,T)+pagemtimes(A(:,2:end,:),'none',A(:,2:end,:),'ctranspose');
    if isnan(w_true)
        w_true = lcmvweights([A(:,1,1), A(:,1,T)], [1; 1], mean(Cx_true,3));
    end
    w_ini = lcmvweights([a_1_ini, a_T_ini], [1; 1], Cx); % LCMP
    B = cat(2,A(2:end,1,:),-A(1,1,:).*eye(d-1));
   
    ini_vals.a_1_ini = a_1_ini;
    ini_vals.a_T_ini = a_T_ini;
    ini_vals.w_ini = w_ini;
   
    true_vals.A = A;
    true_vals.B = B;
    true_vals.a_1_true = A(:,1,1);
    true_vals.a_T_true = A(:,1,end);
    true_vals.w_true = w_true;
    true_vals.sigma_soi = sigma_soi;
    true_vals.Cz = pagemtimes(pagemtimes(B,Cx_true),'none',B,'ctranspose');

    X   = reshape(X, d, N);
    S   = reshape(S, d, N);
    I   = reshape(I, d, N);
    soi = reshape(soi, 1, N);
end
