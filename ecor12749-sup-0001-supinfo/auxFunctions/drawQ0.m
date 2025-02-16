function [empty,S,Q0] = drawQ0(restr,phi,opt)
% Functions attempts to draw Q from the space of orthonormal matrices
% satisfying identifying restrictions.
% Function assumes a particular set of restrictions on the covariances
% between a single proxy variable and the first structural shock.
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form parameters
% - opt: structure containing options

signRestr = restr.signRestr;
vma = phi.vma;
Sigmatr = phi.Sigmatr;
Sigmatrinv = phi.Sigmatrinv;
Sigmatrinvp = Sigmatrinv';

n = size(Sigmatr,1); % Number of endogenous variables in VAR
s = size(signRestr,1); % Number of sign restrictions (excluding normalisations)

%% Construct matrix representing sign restrictions S(phi)*vec(Q) >= 0.
S = zeros(s,n*n);

for ii = 1:s
    
    if signRestr(ii,5) == 1 % Sign restriction on impulse response
        
        S(ii,((signRestr(ii,2)-1)*n+1):signRestr(ii,2)*n) = ...
            vma(signRestr(ii,1),:,signRestr(ii,3)+1)*signRestr(ii,4);
        
    elseif signRestr(ii,5) == 2 % Sign restriction on A0
        
        S(ii,((signRestr(ii,1)-1)*n+1):signRestr(ii,1)*n) = ...
            Sigmatrinv(:,signRestr(ii,2))'*signRestr(ii,4);
        
    end
    
end

%% Add restrictions on covariance between proxy and first shock.
% Restrictions are that (1) proxy is positively correlated with first shock
% and (2) correlation between proxy and first shock is larger in absolute
% value than correlation with other shocks. This is equivalent to:
% 1) d'*Sigma_tr*q_1 >= 0 
% 2) d'*Sigma_tr*(q_1 - q_j) >=0 and d'*Sigma_tr(q_1 + q_j) >= 0 for 1<j<=n

if restr.proxyRestr == 1
    
    d = phi.d'*Sigmatr; 
    S_proxy = zeros(2*n,n*n);   
    S_proxy(:,1:n) = repmat(d,[2*n,1]);
  
    for jj = 2:n
        
        S_proxy(jj,((jj-1)*n+1):jj*n) = -d;
        S_proxy(n+jj,((jj-1)*n+1):jj*n) = d;
        
    end
      
    S = [S; S_proxy]; % Combine with other sign restrictions
      
end

%% Attempt to draw Q satisfying sign restrictions.

iter = 0;
empty = 1;

while iter <= opt.maxDraw && empty == 1

    iter = iter + 1;

    % Draw nxn orthonormal matrix from uniform distribution.
    [Q0,~] = qr(randn(n),0);

    % Normalise diagonal elements of A0 to be positive.
    Q0 = (sign(diag(Sigmatrinvp*Q0))').*Q0;
    
    % Check whether proposed draw satisfies sign restrictions.
    empty = any(S*Q0(:) < 0);

end

end