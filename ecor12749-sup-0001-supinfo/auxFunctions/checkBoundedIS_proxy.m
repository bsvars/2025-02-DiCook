function unbounded = checkBoundedIS_proxy(S,phi,opt)
% Function checks whether the identified set for the normalising variable
% includes zero and thus whether the identified sets for the impulse
% responses to a unit shock are potentially unbounded given a set of 
% traditional and proxy-based sign restrictions.
% Inputs:
% - S: matrix representing coefficients of sign restrictions (including
% proxy-based restrictions)
% - phi: structure containing reduced-form parameters
% - opt: structure containing options

Sigmatrinv = phi.Sigmatrinv;

n = size(Sigmatrinv,1); % Number of endogenous variables in VAR

% Add sign restriction on impact response of normalising variable to first 
% shock to set of zero restrictions.
F = S(1,1:n);
S = S(2:end,:);

%% Attempt to draw Q satisfying restrictions.

iter = 0;
empty = 1;

while iter <= opt.maxDraw && empty == 1
    
    iter = iter + 1;

    % Draw Q from space of orthonormal matrices satisfying zero 
    % restrictions.

    % Generate vector of standard normal random variables.
    z = randn(n,1); 

    Qtilde = zeros(n);
    
    % Compute residual from linear projection of z on F.
    [q,r] = qr(F',0);
    Qtilde(:,1) = z-F'*(r \ (q'*z));

    for jj = 2:n % For each column of Q
       
       % Matrix of coefficients representing orthogonality restrictions.
       regMat = Qtilde(:,1:jj-1);
        
       % Generate vector of independent standard normal random variables.
       z = randn(n,1); 

       % Compute residual from linear projection of z on regMat.
       [q,r] = qr(regMat,0);
       Qtilde(:,jj) = z - regMat*(r \ (q'*z));

    end    
    
    % Normalise diagonal elements of A0 to be positive. Note that Matlab 
    % is implicitly expanding arrays to be compatible with elementwise 
    % array operations.
    Q = ((sign(diag(Sigmatrinv'*Qtilde))').*Qtilde)./vecnorm(Qtilde);
    
    % Check whether proposed draw satisfies sign restrictions.
    empty = any(S*Q(:) < 0);
        
end

unbounded = 1-empty;

end