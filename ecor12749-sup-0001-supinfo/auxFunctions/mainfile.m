%% Construct data for estimating VAR for y_t.
YY = data(opt.p+1:end,:); % y_t 
XX = lags(data,1:opt.p); % Matrix of regressors in VAR for y_t
XX = XX(opt.p+1:end,:); % Drop initial missing observations
if opt.const == 1 % Add constant to matrix of regressors
    XX = [XX ones(size(XX,1),1)];
end
% Add exogenous variables to matrix of regressors.
XX = [XX exog(opt.p+1:end,:)]; 

n = size(YY,2); % Number of endogenous variables
m = size(XX,2); % Number of parameters in each equation for y_t
opt.nExog = opt.const+size(exog,2); % Number coefficients on exogenous variables
T = length(YY); % Number of observations used in estimating VAR

%% Conduct posterior inference on impulse responses.
% Assumes improper (Jeffreys) prior, which generates normal-inverse-Wishart
% posterior.
phiHat.B = (XX'*XX)\XX'*YY;
phiHat.S = (YY - XX*phiHat.B)'*(YY-XX*phiHat.B);
phiHat.Sigma = (1/T)*phiHat.S; % Variance matrix of VAR innovations
phiHat.P = (XX'*XX)\eye(m);
phiHat.cholP = chol(phiHat.P,'lower');

% Storage arrays.
B = zeros(n*m,opt.nonEmpty);
Sigmatr = zeros(n,n,opt.nonEmpty);
etaDraws = zeros(opt.nonEmpty,opt.H+1,length(opt.ivar));
etalb = etaDraws;
etaub = etaDraws;
coeflb = zeros(opt.nonEmpty,length(opt.ivar));
coefub = coeflb;
unbounded = zeros(opt.nonEmpty,1);
phiDraw = 0; % Counter for no. of draws from posterior of phi
draw = 0; % Counter for no. of draws with nonempty identified set

tic

while draw < opt.nonEmpty

    % Sample from normal-inverse-Wishart posterior.
    phiDraw = phiDraw + 1;
    phi.Sigma = iwishrnd(phiHat.S,T-m);
    phi.Sigmatr = chol(phi.Sigma,'lower');
    phi.Sigmatrinv = phi.Sigmatr\eye(n);
    phi.B = phiHat.B(:) + kron(phi.Sigmatr,phiHat.cholP)*randn(m*n,1);
     
    % Generate coefficients in orthogonal reduced-form VMA representation.
    [phi.vma,~] = genVMA(phi,opt);
    
    % Check if identified set is nonempty using approach from
    % Giacomini, Kitagawa and Volpicella (2021).
    [empty,restr.F,restr.S] = checkEmptyIS_GKV(restr,phi);
    
    if empty == 0
    
        draw = draw + 1;
        
        % Store parameters.
        B(:,draw) = phi.B;
        Sigmatr(:,:,draw) = phi.Sigmatr;
        
        % Draw q from uniform distribution over space of unit-length
        % vectors satisfying restrictions.
        q0 = drawq(restr,phi);
        
        % Use draw to compute impulse responses; this represents a draw of
        % the impulse responses given a conditionally uniform prior.
        
        for hh = 1:opt.H+1 % For each horizon
    
            etaDraws(draw,hh,:,:) = phi.vma(opt.ivar,:,hh)*q0;

        end
        
        % Check whether identified sets of impulse responses to a unit
        % shock are potentially unbounded.
        unbounded(draw) = checkBoundedIS_Read(restr.F,restr.S);

        % Compute bounds of identified set for impulse responses using
        % approach from Gafarov, Meier and Montiel-Olea (2018).
        [etalb(draw,:,:),etaub(draw,:,:)] = analyticalBounds(restr,phi,opt);
       
        % Compute bounds of identified set for structural coefficients using
        % modification of above approach.
        [coeflb(draw,:),coefub(draw,:)] = analyticalBounds_coef(restr,phi,opt);
       
        if mod(opt.nonEmpty-draw,opt.dispIter) == 0
              
         fprintf('\n%d draws with non-empty identified set remaining...',...
         opt.nonEmpty-draw);
               
        end
             
    end

end

runTime = toc;

% Compute posterior mean under single prior.
etaMean = permute(mean(etaDraws,1),[2 3 1]);

% Compute highest posterior density intervals under single prior.
[etaHpdlb,etaHpdub] = highestPosteriorDensity(etaDraws,opt);

% Compute set of posterior means.
etaMeanlb = permute(mean(etalb,1),[2 3 1]);
etaMeanub = permute(mean(etaub,1),[2 3 1]);
coefMeanlb = mean(coeflb,1);
coefMeanub = mean(coefub,1);

% Compute robust credible intervals.
[etaCredlb,etaCredub] = credibleRegion(etalb,etaub,opt);

% Compute posterior plausibility of restrictions (posterior probability
% that identified set is nonempty).
postPlaus = draw/phiDraw;