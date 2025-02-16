function q = drawq(restr,phi)
% Function draws q from the space of orthonormal matrices satisfying 
% the identifying restrictions. This function assumes that the 
% restrictions constrain the first column of Q only and that the identified
% set is nonempty.
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form VAR parameters

S = restr.S(1:end,:); % Drop sign normalisation from sign restrictions
F = restr.F;
Sigmatrinv = phi.Sigmatrinv;

n = size(Sigmatrinv,1); % Number of variables in VAR

flag = 0;

while flag == 0

    % Draw nx1 vector of independent standard normal random variables.
    z = randn(n,1);
    
    if isempty(F)
        q = z;
    else
        % Compute residual from linear projection of z on F'. Resulting 
        % vector satisfies zero restrictions.
        [q,r] = qr(F',0);
        q = z-F'*(r \ (q'*z));
    end

    % Normalise vector to satisfy sign normalisation and have unit length.
    q = sign(Sigmatrinv(:,1)'*q)*(q./norm(q));

    % Check if sign restrictions satisfied given draw of q. If so,
    % terminate while loop. Otherwise, increment counter.
    flag = all(S*q >= 0);
        
end

end