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

%% Construct data for estimating 'first-stage' regression.
% Drop initial proxy observations to align with other variables.
proxy = proxy(opt.p+1:end);
% Keep observations where proxy is nonmissing.
proxy = proxy(~isnan(proxy));
ZZ_fs = [YY(~isnan(proxy),:) XX(~isnan(proxy),:)];

m_fs = m+n; % Number of parameters in first-stage regression
T_fs = length(proxy); % Number of observations in first-stage regression

%% Conduct posterior inference on impulse responses.
% Assumes improper (Jeffreys) prior, which generates normal-inverse-Wishart
% posterior.
phiHat.B = (XX'*XX)\XX'*YY;
phiHat.S = (YY - XX*phiHat.B)'*(YY - XX*phiHat.B);
phiHat.Sigma = (1/T)*phiHat.S; % Variance matrix of VAR innovations
phiHat.P = (XX'*XX)\eye(m);
phiHat.cholP = chol(phiHat.P,'lower');

phiHat.D = (ZZ_fs'*ZZ_fs)\ZZ_fs'*proxy;
phiHat.Omega = (proxy - ZZ_fs*phiHat.D)'*(proxy - ZZ_fs*phiHat.D);
phiHat.P_fs = (ZZ_fs'*ZZ_fs)\eye(m_fs);
phiHat.cholP_fs = chol(phiHat.P_fs,'lower');

% Storage arrays.
B = zeros(n*m,opt.nonEmpty);
Sigmatr = zeros(n,n,opt.nonEmpty);
d = zeros(n,opt.nonEmpty);
etaDraw = zeros(opt.nonEmpty,opt.H+1,length(opt.ivar));
etalb = etaDraw;
etaub = etaDraw;
relDraws = etaDraw;
rellb = etaDraw;
relub = etaDraw;
gammalb = zeros(opt.nonEmpty,1);
gammaub = gammalb;
unbounded = zeros(opt.nonEmpty,1);
phiDraw = 0; % Counter for no. of draws from posterior of phi
draw = 0; % Counter for no. of draws with nonempty identified set

tic

while draw < opt.nonEmpty

    % Sample from normal-inverse-Wishart posterior for reduced-form VAR
    % parameters.
    phiDraw = phiDraw + 1;
    phi.Sigma = iwishrnd(phiHat.S,T-m);
    phi.Sigmatr = chol(phi.Sigma,'lower');
    phi.Sigmatrinv = phi.Sigmatr\eye(n);
    phi.B = phiHat.B(:) + kron(phi.Sigmatr,phiHat.cholP)*randn(m*n,1);
     
    % Generate coefficients in orthogonal reduced-form VMA representation.
    [phi.vma,~] = genVMA(phi,opt);
    
    % Sample from normal-inverse-Wishart posterior for first-stage
    % regression parameters.
    phi.Omega = iwishrnd(phiHat.Omega,T_fs-m_fs);
    phi.Omegatr = chol(phi.Omega,'lower');
    phi.D = phiHat.D(:) + kron(phi.Omegatr,phiHat.cholP_fs)*randn(m_fs,1);
    % Extract coefficients on y_t in first-stage regression.
    phi.d = phi.D(1:n);
    
    % Check if identified set is nonempty using rejection sampling.
    [empty,restr.S,Q0] = drawQ0(restr,phi,opt);
    
    if empty == 0
    
        draw = draw + 1;
        
        % Store parameters.
        B(:,draw) = phi.B;
        Sigmatr(:,:,draw) = phi.Sigmatr;
        d(:,draw) = phi.d;
        
        % Use value of Q obtained above to compute impulse responses
        % (i.e. draw from conditionally uniform prior).
        for hh = 1:opt.H+1 % For each horizon

            etaDraw(draw,hh,:) = phi.vma(opt.ivar,:,hh)*Q0(:,1);
            
            % Compute responses to unit shock.
            relDraws(draw,hh,:) = etaDraw(draw,hh,:)./etaDraw(draw,1,1);

        end        

        % Check whether identified sets of impulse responses to a unit
        % shock are unbounded.
        unbounded(draw) = checkBoundedIS_proxy(restr.S,phi,opt);
        
        % Compute bounds of identified set using rejection sampling. 
        out = approximateBounds(restr.S,phi,opt);
        etalb(draw,:,:) = out.rmin;
        etaub(draw,:,:) = out.rmax;
        rellb(draw,:,:) = out.hmin;
        relub(draw,:,:) = out.hmax;
        gammalb(draw) = out.gmin;
        gammaub(draw) = out.gmax;
            
      if mod(opt.nonEmpty-draw,opt.dispIter) == 0
              
       fprintf('\n%d draws with non-empty identified set remaining...',...
       opt.nonEmpty-draw);
       save('proxyRestr_temp.mat');
               
      end
             
    end

end

runTime = toc;

% Compute posterior mean of impulse response under single prior.
etaMean = permute(mean(etaDraw,1),[2 3 1]);
% Compute posterior median of impulse response under single prior.
etaMed = permute(median(etaDraw,1),[2 3 1]);

% Compute highest posterior density intervals under single prior.
[etaHpdlb,etaHpdub] = highestPosteriorDensity(etaDraw,opt);

% Compute sets of posterior means for impulse responses.
etaMeanlb = permute(mean(etalb,1),[2 3 1]);
etaMeanub = permute(mean(etaub,1),[2 3 1]);
% Compute sets of posterior medians for impulse responses.
etaMedlb = permute(median(etalb,1),[2 3 1]);
etaMedub = permute(median(etaub,1),[2 3 1]);

% Robust credible interval for impulse response is smallest volume credible
% region, as in GK.
[etaCredlb,etaCredub] = credibleRegion(etalb,etaub,opt);

% Compute posterior median of impulse response to fixed-size shock.
relMed = permute(median(relDraws,1),[2 3 1]);

% Compute highest posterior density intervals for responses to fixed-size
% shock.
[relHpdlb,relHpdub] = highestPosteriorDensity(relDraws,opt);

% Compute sets of posterior medians for responses to fixed-size shock.
relMedlb = permute(median(rellb,1),[2 3 1]);
relMedub = permute(median(relub,1),[2 3 1]);

% Compute set of posterior means for correlation between proxy and shock.
gammaMeanlb = mean(gammalb)/std(proxy);
gammaMeanub = mean(gammaub)/std(proxy);

% Compute equi-tailed robust credible intervals for covariance.
gammaCredlb = prctile(gammalb,(100-95)/2)/std(proxy);
gammaCredub = prctile(gammaub,100-(100-95)/2)/std(proxy);

% Compute posterior plausibility of restrictions (posterior probability
% that identified set is nonempty).
postPlaus = draw/phiDraw;