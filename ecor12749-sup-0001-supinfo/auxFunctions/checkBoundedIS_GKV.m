function unbounded = checkBoundedIS_GKV(F,S)
% Function determines whether the identified sets for the impulse responses 
% to a unit shock are unbounded by augmenting the set of zero restrictions
% with a binding sign restriction on the impact response of the first
% variable to the first shock and checking whether the associated
% 'identified set' for the impact response of the first variable is
% nonempty.
% The specific approach to checking whether the implied identified set is
% nonempty is from Giacomini, Kitagawa and Volpicella (2021).
% Specifically, consider all possible combinations of active sign 
% restrictions and check whether there exists a vector satisfying these
% constraints and any inactive constraints. 
% This function assumes that the restrictions constrain the first column 
% of Q only and that the first sign restriction restricts the impact
% response of the first variable to the first shock to be nonnegative.
% Inputs:
% - F: matrix containing cofficients of active restrictions
% - S: matrix containing coefficients of nonactive sign restrictions

% Add sign restriction on impact response of first variable to first shock
% to set of zero restrictions,
F = [F; S(1,:)];
S = S(2:end,:);

n = size(F,2); % Number of variables in VAR
f = size(F,1); % Number of zero restrictions (plus additional binding sign restriction)
s = size(S,1); % Number of remaining sign restrictions (including normalisation)

%% Check for existence of vector satisfying restrictions.

% If s + f <= n, there always exists a vector satisfying the restrictions.
if s + f <= n 
    unbounded = 1;
    return
end

% Number of active sign restrictions to consider.
sTilde = n-f-1;
% Generate matrix containing indices representing all possible 
% combinations of sTilde columns of S.
sComb = nchoosek(1:s,sTilde);
% Each row of sComb contains n-f-1 indices representing a combination of
% sign restrictions (indexing the rows of S).
rNum = size(sComb,1); % Number of possible combinations

for rr = 1:rNum % For each possible combination of active constraints

    % Construct matrix containing active constraints.
    Z = [F; S(sComb(rr,:),:)];
    % Construct matrix containing non-active sign restrictions.
    Sna = S(setdiff(1:end,sComb(rr,:)),:);
    % Compute an orthonormal basis for null space of Z, which is a
    % one-dimensional vector since rank(Z) = n-1.
    q = null(Z);
    % Check if q or -q satisfy the non-active sign restrictions.
    unbounded = (all(Sna*q >= 0) || all(-Sna*q >= 0));

    if unbounded == 1 % Conclude identified sets to unit shock are unbounded
        
        return % Terminate function

    end

end
        
end