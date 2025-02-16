% Estimate the effects of monetary policy in Australia using an SVAR with
% sign restrictions on impulse responses, as in Uhlig (2005), plus sign
% restrictions on the systematic component of monetary policy, as in Arias,
% Caldara and Rubio-Ramirez (2019).
% Restriction (5) in paper.

clear variables
close all

oldFolder = pwd;
addpath('auxFunctions');
cd results
figDir = pwd;
cd(oldFolder);

%% Import data.
% The variables are: 
% 1) cash rate (CASH)
% 2) real GDP (GDP)
% 3) trimmed mean CPI (CPI)
% 4) nominal trade-weighted exchange-rate index (TWI) 
% 5) goods and services terms of trade (TOT) 
% 6) US real GDP (USGDP)
% 7) federal funds rate (FFR)
% 8) Monetary policy shock proxy from Beckers (2020) (PROXY)
% 9) Money market spread (MMS)
% 10) Unemployment rate (UR)
DATA = readtable('VARData.xlsx','Sheet','Data');

% Extract date variable and convert into date representation.
date = table2array(DATA(:,1));

% Select sample period (1984Q1:2019Q4).
DATA = table2array(DATA(date>=datetime(1984,3,1) & date<datetime(2020,1,1),2:end));

% Gather endogenous variables and apply desired transformations. Note that
% variable whose shock is of interest (i.e. cash rate) should be first.
data = [DATA(:,1), 100*log(DATA(:,2:4))];

% Set endogenous variable names.
varnames = {'Cash Rate','Real GDP','CPI','TWI'};

%% Options.
opt.p = 4; % No. of lags in VAR
opt.const = 1; % const = 1 if constant in VAR, = 0 otherwise
opt.ivar = 1:4;  % Indices of variables of interest
opt.cumIR = []; % Indices of variables for cumulative impulse responses 
opt.H = 20; % Terminal horizon for impulse responses
opt.nonEmpty = 1000; % No. of draws from posterior of phi with non-empty identified set
opt.aalpha = 68; % Credibility level (%) for credible intervals
opt.dispIter = 100; % Print number of draws remaining every dispIter draws
opt.gridLength = 1000; % Size of grid used when computing credible intervals

% Declare exogenous variables (other than constant; leave empty if none).
exog = [100*log(DATA(:,5:6)), DATA(:,7)];
% Including contemporaneous value and four lags of terms of trade, US GDP
% and federal funds rate as exogenous variables.
exog = lags(exog,0:opt.p);

clear DATA RAW

%% Input identifying restrictions.
% Each row of signRestr contains a vector (i,j,h,s,t) representing a
% 'traditional' sign restriction, where t is the type of restriction:
% t = 1: the impulse response of the ith variable to the jth shock at the 
% hth horizon is nonnegative (s = 1) or nonpositive (s = -1).
% t = 2: the (ij)th element of A0 is nonnegative (s = 1) or nonpositive 
% (s = -1). 
% Functions assume that the first sign restriction restricts the impact
% response of the first variable to the first shock to be nonnegative.
% signRestr = []; % No sign restrictions on impulse responses
restr.signRestr = ...
      [1 1 0 1 1; % Response of IBOCR to monetary policy shock on impact is nonnegative
       1 2 0 -1 2; % Coefficient on real GDP is nonpositive
       1 3 0 -1 2; % Coefficient on CPI is nonpositive
       1 1 1 1 1; % As above after one quarter
       1 1 2 1 1; % As above after two quarters
       1 1 3 1 1; % As above after three quarters
       3 1 0 -1 1; % Response of CPI to monetary policy shock is nonpositive
       3 1 1 -1 1; % As above after one quarter
       3 1 2 -1 1; % As above after two quarters
       3 1 3 -1 1]; % As above after three quarters
  
% Each row of eqRestr contains a vector (i,j,t) representing a 
% particular equality restriction, where t is the type of restriction:
% t = 1: the (ij)th element of A0 is zero
% t = 2: the (ij)th element of A0^(-1) is zero
restr.eqRestr = []; % No zero restrictions
   
%% Conduct (robust) posterior inference.
rng(19061987); % Set seed for random number generator
mainfile;

% Compute posterior probabilities of output falling after two years.
posteriorProb = mean(etaDraws(:,9,2) <0);
fprintf('\nPosterior probability that output falls after two years is %1.4g\n',...
    posteriorProb);
lowerProb = mean(etaub(:,9,2) < 0);
fprintf('\nLower posterior probability that output falls after two years is %1.4g\n',...
    lowerProb);

% Ensure that credible intervals respect sign restrictions on impulse
% responses.
for hh = 1:4
    if etaCredlb(hh,1) < 0
        etaCredlb(hh,1) = 0;
    end
    if etaCredub(hh,3) > 0
        etaCredub(hh,3) = 0;
    end 
end

cd(figDir);

% Create table of results.
for ii = 1:length(opt.ivar)
    
    TT = table((0:opt.H)',etaMean(:,ii),etaHpdlb(:,ii),etaHpdub(:,ii),...
            etaMeanlb(:,ii),etaMeanub(:,ii),etaCredlb(:,ii),etaCredub(:,ii));
    TT.Properties.VariableNames = {'Horizon','Mean','HPDLB','HPDUB',...
        'MeanLB','MeanUB','RCLB','RCUB'};
    writetable(TT,'FigureData.xlsx','Sheet',strcat(varnames{ii},'_signPlusSys'));
    
end

cd(oldFolder);

save('signPlusSys_results.mat');

%% Impose additional restrictions on response to exchange rate.
% Each row of signRestr contains a vector (i,j,h,s,t) representing a
% 'traditional' sign restriction, where t is the type of restriction:
% t = 1: the impulse response of the ith variable to the jth shock at the 
% hth horizon is nonnegative (s = 1) or nonpositive (s = -1).
% t = 2: the (ij)th element of A0 is nonnegative (s = 1) or nonpositive 
% (s = -1). 
restr.signRestr = ...
    [restr.signRestr;
     4 1 0 1 1; % Impact response of TWI is nonnegative
     4 1 1 1 1; % As above after one quarter
     4 1 2 1 1; % As above after two quarters
     4 1 3 1 1; % As above after three quarters
     1 4 0 1 2]; % Coefficient on TWI is nonnegative
 
rng(19061987); % Set seed for random number generator
mainfile; % Conduct (robust) posterior inference

% Compute posterior probabilities of output falling after two years.
posteriorProb = mean(etaDraws(:,9,2) <0);
fprintf('\nPosterior probability that output falls after two years is %1.4g\n',...
    posteriorProb);
lowerProb = mean(etaub(:,9,2) < 0);
fprintf('\nLower posterior probability that output falls after two years is %1.4g\n',...
    lowerProb);

% Ensure that credible intervals respect sign restrictions on impulse
% responses.
for hh = 1:4
    if etaCredlb(hh,1) < 0
        etaCredlb(hh,1) = 0;
    end
    if etaCredub(hh,3) > 0
        etaCredub(hh,3) = 0;
    end
    if etaCredlb(hh,4) < 0
        etaCredlb(hh,4) = 0;
    end    
end

cd(figDir);

% Create table of results.
for ii = 1:length(opt.ivar)
    
    TT = table((0:opt.H)',etaMean(:,ii),etaHpdlb(:,ii),etaHpdub(:,ii),...
            etaMeanlb(:,ii),etaMeanub(:,ii),etaCredlb(:,ii),etaCredub(:,ii));
    TT.Properties.VariableNames = {'Horizon','Mean','HPDLB','HPDUB',...
        'MeanLB','MeanUB','RCLB','RCUB'};
    writetable(TT,'FigureData.xlsx','Sheet',strcat(varnames{ii},'_signPlusSysTWI'));
    
end

cd(oldFolder);

save('signPlusSysTWI_results.mat'); 