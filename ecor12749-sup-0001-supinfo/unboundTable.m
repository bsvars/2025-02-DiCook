% Collate posterior probability that identified set for impact response of
% cash rate includes zero. Also collate posterior probability that
% identified set is nonempty (i.e. posterior plausibility).

clear variables
close all

omega = zeros(6,1);
posteriorPlausibility = zeros(6,1);

load('signRestr_results.mat'); % Load results from Restriction (1)
omega(1) = mean(unbounded);
posteriorPlausibility(1) = postPlaus;

load('signRestrTWI_results.mat'); % Load results from Restriction (2)
omega(2) = mean(unbounded);
posteriorPlausibility(2) = postPlaus;

load('systematicResp_results.mat');
omega(3) = mean(unbounded);
posteriorPlausibility(3) = postPlaus;

load('systematicRespTWI_results.mat');
omega(4) = mean(unbounded);
posteriorPlausibility(4) = postPlaus;

load('signPlusSysTWI_results.mat');
omega(5) = mean(unbounded);
posteriorPlausibility(5) = postPlaus;

load('proxyRestr_results.mat');
omega(6) = mean(unbounded);
posteriorPlausibility(6) = postPlaus;

fprintf('Posterior probability that zero lies within identified set for impact response of cash rate:\n');
omega

fprintf('Posterior plausibility:\n');
posteriorPlausibility