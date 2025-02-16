function empty = chebyCheck(F,S)
% Function checks whether there exists a value of q satisfying the active
% and sign restrictions using Algorithm 1, which extends the approach from
% Amir-Ahmadi and Drautzburg (2020) to allow for zero restrictions.
% Inputs:
% - F: matrix containing cofficients of active restrictions
% - S: matrix containing coefficients of nonactive sign restrictions

n = size(F,2); % Number of variables in VAR
r = size(F,1); % Number of active restrictions
s = size(S,1); % Number of nonactive sign restrictions 

% Options for linear programming routine.
lpoptimOptions = optimoptions('linprog');
lpoptimOptions.Display = 'off';

%% Transform restrictions so that algorithm from AD21 can be applied.
% Compute change of basis matrix.
K = [null(F), null(null(F)')];

% Apply change of basis to coefficient vectors of active restrictions.
Stilde = (K'*S')'; % Note that K^{-1} = K' since K is orthonormal

% Project sign restrictions onto hyperplane generated by zero restrictions.
b = [zeros(r,n-r), eye(r)]';
B = eye(n) - b*((b'*b)\b');
Sbar = (B*Stilde')';

Sbar = Sbar(:,1:(n-r)); % Drop last f elements (which are zero)

%% Solve for Chebyshev centre.
% Find the centre and radius of the largest ball that can be inscribed
% within the intersection of the half-spaces generated by the (projected)
% sign restrictions and the unit cube in n-r dimensions.

A = zeros(s+1,n-r+1);
A(1:end-1,2:end) = -Sbar;

for ii = 1:(s+1)

    A(ii,1) = norm(A(ii,2:end));

end

A(end,1) = -1; % Add constraint that radius is positive

% Add additional restrictions that ball lies inside unit (n-f)-cube.
A = [A; [ones(n-r,1), eye(n-r)]; [ones(n-r,1), -eye(n-r)]];

% Problem is to maximise a*x s.t. A*x <= b, with x = (R,c')'.
a = [-1, zeros(1,n-r)]; % 
b = [zeros(s+1,1); ones(2*(n-r),1)];
x = linprog(a,A,b,[],[],[],[],lpoptimOptions);

% Check if radius is positive (i.e. if set is nonempty).
if x(1) > 0
    
    empty = 0;
    
else
    
    empty = 1;

end

end