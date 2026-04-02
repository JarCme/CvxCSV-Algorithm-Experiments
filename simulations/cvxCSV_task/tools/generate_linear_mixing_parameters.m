function [A, w] = generate_linear_mixing_parameters(d, T, domain, type, angle)

    lambda = linspace(0,1,T);

    switch domain
        case 'complex'
            d_randn = @(dims) crandn(dims);
        case 'real'
            d_randn = @(dims) randn(dims);
        otherwise
            error('Unknown domain.');
    end

    switch type
        case 'model'
            h = d_randn([d-1,1]);
            beta = d_randn(1);  
    
            g_1 = d_randn([d-1,1]);
            g_T = d_randn([d-1,1]);
    
            w = [beta;h];
                
            a = zeros(d,1,T);
            Q = zeros(d,d-1,T);
            g = (1-lambda).*g_1 + (lambda).*g_T;
    
            for t = 1:T
                if T == 1
                    g_t = g_1; 
                else
                    g_t = g(:,t);
                end
                gamma_t = (1-h'*g_t)/conj(beta);
                a_t = [gamma_t;g_t];
    
                Q_t = [h'; (1/gamma_t)*(g_t*h'-eye(d-1))];
                a(:,:,t) = a_t;
                Q(:,:,t) = Q_t;
            end
    
            A = cat(2,a,Q);
            
        case 'model_eye'
            a = eye(d,1);        
            w  = eye(d,1);
            Q  = [zeros(1,d-1);-eye(d-1)];
            A  = zeros(d,d,T);
            A(:,1,:) = repmat(a,1,1,T);
            A(:,2:end,:) = repmat(Q,1,1,T);

        case {'angle_Q_static','angle'}
            a_1 = d_randn([d,1]);
            a_T = randvecangle(a_1,angle);
    
            a = (1-lambda).*a_1 + (lambda).*(a_T); 
            Q = d_randn([d,d-1]); % static over blocks  
            A = cat(2,permute(a,[1,3,2]),repmat(Q,1,1,T));
            w = nan; % true w is not known

        case {'angle_Q_gamma_static','angle_static_gamma'}
            a_1 = d_randn([d,1]);
            a_T = randvecangle(a_1,angle);
            gamma = 1;%crandn(1);
            a_1 = (a_1/a_1(1))*gamma;
            a_T = (a_T/a_T(1))*gamma;
            
    
            a = (1-lambda).*a_1 + (lambda).*(a_T); 
            Q = d_randn([d,d-1]); % static over blocks  
            A = cat(2,permute(a,[1,3,2]),repmat(Q,1,1,T));
            w = nan; % true w is not known

        case 'angle_Q_dynamic'
            a_1 = d_randn([d,1]);
            a_T = randvecangle(a_1,angle);
    
            a = (1-lambda).*a_1 + (lambda).*(a_T); 
            Q = d_randn([d,d-1,T]); % static over blocks  
            A = cat(2,permute(a,[1,3,2]),Q);
            w = nan; % true w is not known
        case 'random' 
            a_1 = d_randn([d,1]);
            a_T = d_randn([d,1]);
            
            a = (1-lambda).*a_1 + (lambda).*(a_T); 
            Q = d_randn([d,d-1]); % static over blocks  
            A = cat(2,permute(a,[1,3,2]),repmat(Q,1,1,T));
            w = nan; % true w is not known
        otherwise 
            error('Unknown type')
    end

end

