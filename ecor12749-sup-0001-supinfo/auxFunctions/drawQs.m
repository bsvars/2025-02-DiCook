function Q = drawQs(S,phi,opt)
% Function attempts to draw Q many times from the space of orthonormal 
% matrices satisfying identifying restrictions.
% Inputs:
% - S: matrix representing coefficients of sign restrictions
% - phi: structure containing reduced-form parameters
% - opt: structure containing options

Sigmatrinvp = phi.Sigmatrinv';

n = size(Sigmatrinvp,1); % Number of endogenous variables in VAR

%% Attempt to draw Q satisfying sign restrictions.

Q = zeros(n,n,opt.Qdraws);

parfor kk = 1:opt.Qdraws
    
    flag = 0;

    while flag == 0

        % Draw nxn orthonormal matrix from uniform distribution.
        [Q0,~] = qr(randn(n),0);
    
        % Normalise diagonal elements of A0 to be positive.
        Q0 = (sign(diag(Sigmatrinvp*Q0))').*Q0;

        % Check whether proposed draw satisfies sign restrictions.
        flag = all(S*Q0(:) >= 0);
        
        if flag == 1
            
            Q(:,:,kk) = Q0;

        end

    end

end

end