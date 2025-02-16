% Estimate the effects of monetary policy in Australia using an SVAR with
% sign restrictions on the covariances between proxies and shocks,
% alongside restrictions on impulse responses (as in Uhlig (2005)) and the
% systematic component of monetary policy, as in Arias, Caldara and 
% Rubio-Ramirez (2019).
% This file replaces real GDP with the unemployment rate.

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
data = [DATA(:,1), DATA(:,10), 100*log(DATA(:,3:4))];

% Extract proxy variable.
proxy = DATA(:,8);

% Set variable names.
varnames = {'Cash Rate','Unemployment','CPI','TWI'};

%% Options.
opt.p = 4; % No. of lags in VAR
opt.const = 1; % const = 1 if constant in VAR, = 0 otherwise
opt.ivar = 1:4;  % Indices of variables of interest
opt.jshock = 1; % Index of structural shock of interest
opt.cumIR = []; % Indices of variables for cumulative impulse responses 
opt.H = 20; % Terminal horizon for impulse responses
opt.nonEmpty = 1000; % No. of draws from posterior of phi with non-empty identified set
opt.aalpha = 68; % Credibility level (%) for credible intervals
opt.dispIter = 100; % Print number of draws remaining every dispIter draws
opt.gridLength = 1000; % Size of grid used when computing credible intervals
opt.maxDraw = 100000; % Maximum number of draws from O(n) when checking emptiness
opt.Qdraws = 5000; % Number of draws of Q when approximating bounds

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
% signRestr = []; % No sign restrictions on impulse responses
restr.signRestr = ...
      [1 2 0 1 2; % Coefficient on unemployment is nonnegative
       1 3 0 -1 2; % Coefficient on CPI is nonpositive
       1 1 0 1 1; % Response of IBOCR to monetary policy shock on impact is nonnegative
       1 1 1 1 1; % As above after one quarter
       1 1 2 1 1; % As above after two quarters
       1 1 3 1 1; % As above after three quarters
       3 1 0 -1 1; % Response of CPI to monetary policy shock is nonpositive
       3 1 1 -1 1; % As above after one quarter
       3 1 2 -1 1; % As above after two quarters
       3 1 3 -1 1]; % As above after three quarters
   
restr.signRestr = ...
    [restr.signRestr;
     4 1 0 1 1; % Impact response of TWI is nonnegative
     4 1 1 1 1; % As above after one quarter
     4 1 2 1 1; % As above after two quarters
     4 1 3 1 1; % As above after three quarters
     1 4 0 1 2]; % Coefficient on TWI is nonnegative   
   
% Each row of eqRestr contains a vector (i,j,t) representing a 
% particular equality restriction, where t is the type of restriction:
% t = 1: the (ij)th element of A0 is zero
% t = 2: the (ij)th element of A0^(-1) is zero
restr.eqRestr = []; % No zero restrictions

% if proxyRestr = 1, impose restrictions on correlation between proxies 
% and shocks. Restrictions are 1) proxy is positively correlated with first
% shock and 2) correlation between proxy and first shock is larger in
% absolute value than correlation wih other shocks 
restr.proxyRestr = 1; 
   
%% Conduct (robust) posterior inference.
rng(19061987); % Set seed for random number generator

mainfile_proxy;

% Compute posterior probabilities of output falling after two years.
posteriorProb = mean(etaDraw(:,9,2) > 0);
fprintf('\nPosterior probability that unemployment rises after two years is %1.4g\n',...
    posteriorProb);
lowerProb = mean(etalb(:,9,2) > 0);
fprintf('\nLower posterior probability that unemployment rises after two years is %1.4g\n',...
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
    
    if relHpdub(hh,3) > 0
        relHpdub(hh,3) = 0;
    end
    if relHpdlb(hh,4) < 0
        relHpdlb(hh,4) = 0;
    end    
end
    
cd(figDir);

% Create table of results for impulse responses.
for ii = 1:length(opt.ivar)
    
    TT1 = table((0:opt.H)',etaMean(:,ii),etaHpdlb(:,ii),etaHpdub(:,ii),...
            etaMeanlb(:,ii),etaMeanub(:,ii),etaCredlb(:,ii),etaCredub(:,ii));
    TT1.Properties.VariableNames = {'Horizon','Mean','HPDLB','HPDUB',...
        'MeanLB','MeanUB','RCLB','RCUB'};
    writetable(TT1,'FigureData.xlsx','Sheet',strcat(varnames{ii},'_proxyRestr_UR'));
    
end

% Create table of results for relative impulse responses.
for ii = 1:length(opt.ivar)
    
    TT2 = table((0:opt.H)',relMed(:,ii),relHpdlb(:,ii),relHpdub(:,ii),...
        relMedlb(:,ii),relMedub(:,ii));
    TT2.Properties.VariableNames = ...
        {'Horizon','Median','HPDLB','HPDUB','MedLB','MedUB'};
    writetable(TT2,'FigureData.xlsx','Sheet',strcat(varnames{ii},'_proxyRestrRel_UR'));
    
end

cd(oldFolder);

% Compute and tabulate posterior lower probability of different events.
cd(figDir);

hors = [0 4 8 12 16]' + 1; % Horizons of interest 
thresh = -(0.25:0.25:1); % Thresholds of interest
threshUR = -thresh; 

lowerProbs = zeros(length(hors),length(thresh),2);

for ii = 1:length(hors)
        
        lowerProbs(ii,:,1) = mean(rellb(:,hors(ii),2) >= threshUR);
        lowerProbs(ii,:,2) = mean(relub(:,hors(ii),3) <= thresh);
        
end

TT4 = table(hors-1,100*lowerProbs(:,:,1),100*lowerProbs(:,:,2));
TT4.Properties.VariableNames = {'Horizon','Unemployment','Prices'};
writetable(TT4,'TableData.xlsx','Sheet','LowerProbs_UR');

upperProbs = zeros(length(hors),length(thresh),2);

for ii = 1:length(hors)
        
        upperProbs(ii,:,1) = mean(relub(:,hors(ii),2) >= threshUR);
        upperProbs(ii,:,2) = mean(rellb(:,hors(ii),3) <= thresh);
        
end

TT5 = table(hors-1,100*upperProbs(:,:,1),100*upperProbs(:,:,2));
TT5.Properties.VariableNames = {'Horizon','Unemployment','Prices'};
writetable(TT5,'TableData.xlsx','Sheet','UpperProbs_UR');

cd(oldFolder);

save('proxyRestr_UR_results.mat');