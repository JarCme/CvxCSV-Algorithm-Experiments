function [w, a, shat, NumIt,w_all] = bogive_w(x, wini, Nb, mu, nonln, maxit, precond, pdf_alpha, pdf_gamma)
% Block-Wise Independent Vector/Component Extraction by Orthogonally-Constrained Gradient
% Algorithm - Gaussian Background, constant separating vector (w) over blocks 

epsilon = 0.0000001; % stopping threshold

[d, N, K] = size(x);
if isreal(x), realvalued = 1; else, realvalued = 0; end

if nargin < 7
    precond = repmat(eye(d),[1 1 K]);
end

if nargin < 6
    maxit = 1000;
end

if nargin < 5
    nonln = 'rati';
end

if nargin < 4
    mu = 0.1;
end

if nargin < 3
    Nb = N; % conventional static ICE/IVE
end

if nargin < 2
    if realvalued
        wini = randn(d, K);
    else
        wini = crandn(d, K);
    end
end


w_all = zeros(d,1,K,maxit+1);

%%%%%%%%%%%% Blocks
x = x - mean(x,2);

M = max(floor(N/Nb),1); % the number of blocks
Nb = floor(N/M);
N = Nb*M;
x = x(:,1:N,:); % the incomplete block is neglected

X = permute(reshape(permute(x,[2 1 3]), Nb, M, d, K),[3 1 4 2]);
Cb_in = mtimesx(X,X,'C')/Nb;
C_in = mean(Cb_in,4); % the covariance of the whole data

%%%%%%%%%%%% Preconditioning (preprocessing)
if isempty(precond)
    precond = zeros(d,d,K);
    for k = 1:K
        precond(:,:,k) = sqrtm(inv(C_in(:,:,k)));
    end
end
x = mtimesx(precond, x);
wini = permute(wini,[1 3 2]);
for k = 1:K
    wini(:,1,k) = (precond(:,:,k)')\wini(:,1,k);
end

X = permute(reshape(permute(x,[2 1 3]), Nb, M, d, K),[3 1 4 2]);

Cb = mtimesx(X,X,'C')/Nb;
C = mean(Cb,4);

% OG initialization & scaling
w = wini;
% for k = 1:K % MPDR steered in the direction given by aini(:,1,k)
%     aux = C(:,:,k)\aini(:,1,k);
%     w(:,1,k) = aux/(aini(:,1,k)'*aux);
% end
w_all(:,:,:,1) = w;
w = w./sqrt(mtimesx(mtimesx(w,'C',C),w)); % output scale normalized to one
%aux = mtimesx(Cb,w);
%a = aux./mtimesx(w,'C',aux); % mixing vectors re-computed

NumIt = 0;
crit = 0;
%while crit < 1-epsilon && NumIt < maxit
while NumIt < maxit

    NumIt = NumIt + 1;
    
   % sigma2 = sum(conj(w).*mtimesx(Cb,w),1); % variance of SOI on blocks
    
    a = mtimesx(Cb,w); 
    sigma2 = sum(conj(w).*a,1); % variance of SOI on blocks
    a = a./sigma2; % mixing vectors obtained through the orthogonal constraint
    wold = w;
    
    soi = mtimesx(w,'C',X);
    sigma = sqrt(sigma2); 
    soin = soi./sigma; % block-normalized SOI 
    if realvalued
        [psi, ~] = realnonln(soin, nonln);
    else
        switch lower(nonln)
            case 'cggd'
                [psi, ~] = complexnonln(soin, nonln, pdf_alpha, pdf_gamma);
            otherwise
                [psi, ~] = complexnonln(soin, nonln);
        end
    end
    
    xpsi = (mtimesx(X,psi,'T')/Nb)./sigma;
    nu = mtimesx(w,'c',xpsi);
    gradw = a - xpsi./nu;
    gradw = sum(gradw,4); % gradient
    w = w + mu*gradw;
    w_all(:,:,:,NumIt+1) = w;
    w = w./sqrt(sum(mtimesx(C,w).*conj(w),1));
    crit = min(squeeze(abs(mtimesx(w,'C',wold))./...
        sqrt(sum(w.*conj(w),1))./sqrt(sum(wold.*conj(wold),1))));
end


% least-squares rescaling by projecting on the first input channel
aux = conj(mtimesx(C_in(1,:,:),mtimesx(precond,'C',w))./mtimesx(w,'C',mtimesx(C,w)));
shat = permute(mtimesx(aux.*w,'C',x),[3 2 1]); 
w = aux.*mtimesx(precond,'C',w);
aux = mtimesx(Cb_in,w);
a = permute(aux./mtimesx(w,'C',aux),[1 4 3 2]);
%w = permute(w,[1 3 2]);
end

%%%%%%%%%%% helping functions

function [psi, psipsi] = realnonln(s,nonln)
    if strcmp(nonln,'sign')
        if size(s,3)==1, error('QuickIVE: Nonlinearity "sign" cannot be used in the real-valued case!'); end
        aux = 1./sqrt(sum(s.^2,3));
        psi = s.*aux;
        psipsi = aux.*(1-psi.^2);
    elseif strcmp(nonln,'tanh')
        if size(s,1) > 1
            aux = 1./sqrt(sum(s.^2,3));
            th = tanh(s);
            psi = th.*aux;
            psipsi = aux.*(1 - th.^2 - psi.*aux);
        else
            psi = tanh(s);
            psipsi = 1 - psi.^2;
        end
    elseif strcmp(nonln,'rati')
        aux = 1./(1+sum(s.^2,3));
        psi = s.*aux;
        psipsi = aux - 2*psi.^2;
    elseif strcmp(nonln,'cggd')
        error('not implemented for real case');
    end
end

function [psi, psipsi] = complexnonln(s,nonln,pdf_alpha,pdf_gamma)
    if strcmp(nonln,'sign')
        sp2 = s.*conj(s);
        aux = 1./sqrt(sum(sp2,3));
        psi = conj(s).*aux;
        psipsi = aux.*(1-psi.*conj(psi)/2);
    elseif strcmp(nonln,'rati')
        sp2 = s.*conj(s);
        aux = 1./(1+sum(sp2,3));
        psi = conj(s).*aux;
        psipsi = aux - psi.*conj(psi);     
    elseif strcmp(nonln,'cggd')
        psi = g_pdf(s, pdf_alpha, pdf_gamma);
        psipsi = -1;
    end
end