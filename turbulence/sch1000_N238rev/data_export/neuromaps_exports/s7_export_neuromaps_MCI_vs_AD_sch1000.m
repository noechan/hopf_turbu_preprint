%% Compute MCI–AD turbulence difference (lambda = 0.01, Schaefer-1000)
% Input: full N145 lambda=0.01 table; MCI and AD are selected by Group.
% Output variable: diff_mci_ad_lambda_0_01 (1 x 1000 double)

clear; clc;
sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

% --- Paths ---
in_xlsx = fullfile(P.harmonization_results, ...
    'Turbu_ComBat_ADNI3_HC_MCI_AD_ABeta_lambda_0_01_sch1000_N145.xlsx');
out_mat = fullfile(P.neuromaps_annotations, ...
    'brain_nodes_diff_mci_ad_lambda_0_01_sch1000_ComBat_N145.mat');

% --- Read table ---
T = readtable(in_xlsx);

% Expect columns: PTID, Group, Schaefer_1 ... Schaefer_1000
% Check column names (optional sanity check)
disp(T.Properties.VariableNames(1:10));

% Identify Schaefer columns
isSchaefer = startsWith(T.Properties.VariableNames, 'Schaefer_');
schaeferCols = T.Properties.VariableNames(isSchaefer);

% Convert to numeric array (subjects x parcels)
X = table2array(T(:, schaeferCols));   % size: [Nsubj x 1000]

% --- Define groups ---
% MCI group:   Group == 'MCI_ABpos'
% AD group:   Group == 'AD_ABpos'
mci_idx = strcmp(T.Group, 'MCI_ABpos');
ad_idx = strcmp(T.Group, 'AD_ABpos');

X_MCI = X(mci_idx, :);   % MCI subjects x 1000
X_AD = X(ad_idx, :);   % AD subjects x 1000

fprintf('N(MCI) = %d\n', sum(mci_idx));
fprintf('N(AD) = %d\n', sum(ad_idx));

% --- Compute parcel-wise group means ---
MCI_mean = mean(X_MCI, 1, 'omitnan');   % 1 x 1000
AD_mean = mean(X_AD, 1, 'omitnan');   % 1 x 1000

% --- Difference: MCI - AD ---
diff_mci_ad_lambda_0_01 = MCI_mean - AD_mean;  % 1 x 1000

% Optional: check a few values
disp('First 5 values of diff_mci_ad_lambda_0_01:');
disp(diff_mci_ad_lambda_0_01(1:5));

% --- Save to .mat ---
save(out_mat, 'diff_mci_ad_lambda_0_01');

fprintf('Saved %s with variable diff_mci_ad_lambda_0_01 (size: %d x %d)\n', ...
    out_mat, size(diff_mci_ad_lambda_0_01,1), size(diff_mci_ad_lambda_0_01,2));
