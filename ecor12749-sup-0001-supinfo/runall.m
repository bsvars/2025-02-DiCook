% This script calls the scripts used to generate the main results in the
% paper. Note that it may take a long time to finish running!

main_signRestr
fprintf('Results obtained under Restrictions (1) and (2).\n');
fprintf(sprintf('This took %1.2g hours\n',runTime/(60*60)));

main_systematicResp
fprintf('Results obtained under Restrictions (3) and (4).\n');
fprintf(sprintf('This took %1.2g hours\n',runTime/(60*60)));

main_signPlusSys
fprintf('Results obtained under Restrictions (5).\n');
fprintf(sprintf('This took %1.2g hours\n',runTime/(60*60)));

main_proxyRestr
fprintf('Results obtained under Restriction (6).\n');
fprintf(sprintf('This took %1.2g hours\n',runTime/(60*60)));

main_proxyRestr_FD
fprintf('Results obtained under Restriction (6) with model in first differences.\n');
fprintf(sprintf('This took %1.2g hours\n',runTime/(60*60)));

main_proxyRestr_IT
fprintf('Results obtained under Restriction (6) and inflation-targeting sample.\n');
fprintf(sprintf('This took %1.2g hours\n',runTime/(60*60)));

main_proxyRestr_Cred
fprintf('Results obtained under Restriction (6) with restrictions related to credit spreads.\n');
fprintf(sprintf('This took %1.2g hours\n',runTime/(60*60)));

main_proxyRestr_UR
fprintf('Results obtained under Restriction (6) replacing output with unemployment.\n');
fprintf(sprintf('This took %1.2g hours\n',runTime/(60*60)));

% Generate results in Table 1 (Informativeness of identifying restrictions)
informativenessTable

% Collate results in Table 2 (Posterior probability that identified set
% for impact response of cash rate includes zero)
unboundTable
