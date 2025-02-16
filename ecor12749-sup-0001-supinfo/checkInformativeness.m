% This file uses the reduced-form parameter draws obtained under
% Restriction (6) to impose Restriction (5). The aim is to ascertain
% whether the apparent informativeness of Restriction (6) is being driven
% by the different sets of reduced-form parameters with nonempty identified
% set.

clear variables
close all
clc

oldFolder = pwd;
addpath('auxFunctions');
cd ../Figures 
figDir = pwd;
cd(oldFolder);

load('proxyRestr_results.mat');

etalb_5 = zeros(opt.nonEmpty,opt.H+1,length(opt.ivar));
etaub_5 = etalb_5;

for draw = 1:opt.nonEmpty % For each draw of reduced-form parameters
    
    % Extract relevant reduced-form parameters.
    phi.B = B(:,draw);
    phi.Sigmatr = Sigmatr(:,:,draw);
    % Re-compute reduced-form impulse responses.
    [phi.vma,~] = genVMA(phi,opt);
    
    % Re-compute restriction coefficient matrices (note that Restriction 
    % (5) is the same as Restriction (6) but omits the proxy-based 
    % restrictions).
    [empty,restr.F,restr.S] = checkEmptyIS_GKV(restr,phi);
    
    % Compute bounds of identified set under Restriction (5). 
    [etalb_5(draw,:,:),etaub_5(draw,:,:)] = ...
        analyticalBounds(restr,phi,opt);

end

% Compute set of posterior means.
etaMeanlb_5 = permute(mean(etalb_5,1),[2 3 1]);
etaMeanub_5 = permute(mean(etaub_5,1),[2 3 1]);

% Compute robust credible intervals.
[etaCredlb_5,etaCredub_5] = credibleRegion(etalb_5,etaub_5,opt);

[etaMeanlb(:,2) etaMeanub(:,2) etaMeanlb_5(:,2) etaMeanub_5(:,2)]




