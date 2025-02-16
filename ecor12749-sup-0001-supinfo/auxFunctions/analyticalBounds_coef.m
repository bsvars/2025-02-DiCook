function [rmin,rmax] = analyticalBounds_coef(restr,phi,opt)
% Function uses analytical results from Gafarov, Meier and Montiel-Olea
% (2018) (GMM18) to find the upper and lower bounds of the identified set
% for the coefficients in the structural equation (i.e., the first row of
% A_0 or B^(-1)).
% In notation of GMM18, model is 
% Y_t = c + A_1*Y_t-1 + ... + A_p*Y_t-p + B*eps_t, eps_t ~ N(0,I).
% Restrictions are on columns of B rather than on Q.
% B = A0^(-1) = Sigmatr*Q.
% Inputs:
% - restr: structure representing restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options

F = restr.F;
S = restr.S;
Sigmatr = phi.Sigmatr;
Sigmatrinv = phi.Sigmatrinv;
ivar = opt.ivar;

n = size(Sigmatr,1); % Number of variables in VAR
f = size(F,1); % Number of zero restrictions
s = size(S,1); % Number of sign restrictions (including normalisation)

Sigmatrprime = Sigmatr';
Sigma = Sigmatr*Sigmatrprime;
Sigmainv = Sigma\eye(n);

% Convert restrictions to be consistent with notation from GMM18.
F = (F*Sigmatrinv)';
S = (S*Sigmatrinv)';

% Compute number of combinations of active sign restrictions.
kmax = min(s,n-f-1); % Maximum number of active sign restrictions
rNum = 0;

for kk = 0:kmax

    rNum = nchoosek(s,kk) + rNum;

end

% Storage
fMax = zeros(rNum,length(ivar));
fMin = zeros(rNum,length(ivar));

cbar = 10^8; % Set penalty term

%% No active sign restrictions (k=0).
Z = F; % Active restrictions (equality restrictions only)

M = eye(n)-Sigmatrprime*Z*((Z'*Sigma*Z)\Z'*Sigmatr);

for ii = 1:length(ivar) % For each structural coefficient
    
    CC = Sigmainv(ii,:);

    % Compute value function (up to sign) given active 
    % restrictions
    v = abs(sqrt(CC*Sigmatr*M*Sigmatrprime*CC'));

    if v ~= 0 % If value function non-zero

        % Potential solutions.
        xPlus = Sigmatr*M*Sigmatrprime*CC'./v;
        xMinus = -xPlus;

        fMaxPlus = v - 2*(1-all(S'*xPlus >= 0))*cbar;
        fMaxMinus = -v - 2*(1-all(S'*xMinus >= 0))*cbar;

        fMinPlus = v + 2*(1-all(S'*xPlus >= 0))*cbar;
        fMinMinus = -v + 2*(1-all(S'*xMinus >= 0))*cbar;

        fMax(1,ii) = max(fMaxPlus,fMaxMinus);
        fMin(1,ii) = min(fMinPlus,fMinMinus);          

    elseif v == 0

        % Check whether there exists a point xstar!=0 satisfying the
        % active restrictions in Z and the nonactive sign restrictions.

        if size(Z,2) < n-1

            % Use Chebyshev criterion.
            empty = chebyCheck(Z'*Sigmatr,S'*Sigmatr);

        elseif size(Z,2) == n-1

            xstar = null(Z'*Sigmatr);
            empty = ~(all(S'*Sigmatr*xstar >= 0) | ...
                all(-S'*Sigmatr*xstar >= 0));

        end
        
        if empty == 0

            fMax(1,ii) = 0;
            fMin(1,ii) = 0;

        elseif empty == 1

            fMax(1,ii) = -2*cbar;
            fMin(1,ii) = 2*cbar;

        end

    end

end




%% 0 < k <= kmax (for different numbers of active sign restrictions)

rCount = 1; % Counter for combinations of active restrictions

for kk = 1:kmax

    % Generate vector containing column indices representing all 
    % possible combinations of k of the columns of S.
    sComb = nchoosek(1:s,kk);

    for rr = 1:size(sComb,1) % For each set of active constraints

        rCount = rCount + 1;
        
        % Active restrictions.
        Z = [F S(:,sComb(rr,:))];
        % Non-active sign restrictions.
        Sna = S(:,setdiff(1:end,sComb(rr,:)));
        
        M = eye(n) - Sigmatrprime*Z*((Z'*Sigma*Z)\Z'*Sigmatr);

        for ii = 1:length(ivar) % For each variable of interest

            % Reduced-form IR of ith variable at hth horizon.
            CC = Sigmainv(ii,:);

            % Compute value function (up to sign) given active 
            % restrictions
            v = abs(sqrt(CC*Sigmatr*M*Sigmatrprime*CC')); 

            if v~=0 % If value function non-zero

                % Potential solutions.
                xPlus = Sigmatr*M*Sigmatrprime*CC'./v;
                xMinus = -xPlus;

                fMaxPlus = v - 2*(1-all(Sna'*xPlus >= 0))*cbar;
                fMaxMinus = -v - 2*(1-all(Sna'*xMinus >= 0))*cbar;

                fMinPlus = v + 2*(1-all(Sna'*xPlus >= 0))*cbar;
                fMinMinus = -v + 2*(1-all(Sna'*xMinus >= 0))*cbar;

                fMax(rCount,ii) = max(fMaxPlus,fMaxMinus);
                fMin(rCount,ii) = min(fMinPlus,fMinMinus);

            elseif v == 0

            % Check whether there exists a point xstar!=0 satisfying the
            % active restrictions in Z and the nonactive sign restrictions.

            if size(Z,2) < n-1

                % Use Chebyshev criterion.
                empty = chebyCheck(Z'*Sigmatr,Sna'*Sigmatr);

            elseif size(Z,2) == n-1

                xstar = null(Z'*Sigmatr);
                empty = ~(all(Sna'*Sigmatr*xstar >= 0) | ...
                    all(-Sna'*Sigmatr*xstar >= 0));

            end

            if empty == 0

                fMax(1,ii) = 0;
                fMin(1,ii) = 0;

            elseif empty == 1

                fMax(1,ii) = -2*cbar;
                fMin(1,ii) = 2*cbar;

            end                 

            end

        end

    end

end

% Compute bounds as max/min over fMax/fMin, respectively, for each
% horizon and variable.
rmax = max(fMax,[],1);
rmin = min(fMin,[],1);

end