% Compute informativeness of identifying restrictions following approach in
% Giacomini and Kitagawa (2021).

clear variables
close all

addpath('auxFunctions');

load('signRestr_results.mat'); % Load results from Restriction (1)

rng(19061987);

width = zeros([size(etaMeanlb),6]);

% Compute width of set of posterior means for each impulse response.
width(:,:,1) = etaMeanub-etaMeanlb;

% Replace identifying restrictions with a minimal/nested set of
% restrictions (i.e. sign normaliastion + sign restrictions on impact
% response of cash rate to a monetary policy shock).
restr.signRestr = restr.signRestr(1,:);

% Compute robust Bayes outputs.
mainfile;

width0 = etaMeanub-etaMeanlb;

% Compute width of set of posterior means for each impulse response under
% different sets of restrictions.

load('signRestrTWI_results.mat'); % Load results from Restriction (2)
width(:,:,2) = etaMeanub-etaMeanlb;

load('systematicResp_results.mat');
width(:,:,3) = etaMeanub-etaMeanlb;

load('systematicRespTWI_results.mat');
width(:,:,4) = etaMeanub-etaMeanlb;

load('signPlusSysTWI_results.mat');
width(:,:,5) = etaMeanub-etaMeanlb;

load('proxyRestr_results.mat');
width(:,:,6) = etaMeanub-etaMeanlb;

% Compute informativeness of restrictions (relative to Model 0).
restrInfo = zeros(size(width));

for ii = 1:6

    restrInfo(:,:,ii) = 1 - width(:,:,ii)./width0;

end

% Informativeness (output response at two-year horizons)
squeeze(restrInfo(8,2,:)) 