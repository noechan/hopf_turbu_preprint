clear; clc;

% =========================================================
% Paths
% =========================================================
data_dir = '/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/plots/';
out_csv  = fullfile(data_dir, 'ADNI3_StageBased_Hopf_BoxplotData.csv');

% =========================================================
% File list and variable names
% =========================================================
groups = {'HC_ABneg', 'HC_ABpos', 'MCI_ABpos', 'AD_ABpos'};

all_group = {};
all_infocap = [];
all_suscep = [];

% =========================================================
% Load each stage file and append to long-format table
% =========================================================
for i = 1:numel(groups)

    grp = groups{i};
    mat_file = fullfile(data_dir, ['perturbation_' grp '.mat']);

    S = load(mat_file);

    infocap_var = [grp '_infocap'];
    suscep_var  = [grp '_suscep'];

    if ~isfield(S, infocap_var)
        error('Variable %s not found in %s', infocap_var, mat_file);
    end

    if ~isfield(S, suscep_var)
        error('Variable %s not found in %s', suscep_var, mat_file);
    end

    infocap_vals = S.(infocap_var);
    suscep_vals  = S.(suscep_var);

    % Force column vectors
    infocap_vals = infocap_vals(:);
    suscep_vals  = suscep_vals(:);

    if numel(infocap_vals) ~= numel(suscep_vals)
        error('Different lengths for info capacity and susceptibility in %s', grp);
    end

    n = numel(infocap_vals);

    all_group   = [all_group; repmat({grp}, n, 1)];
    all_infocap = [all_infocap; infocap_vals];
    all_suscep  = [all_suscep; suscep_vals];
end

% =========================================================
% Create table
% =========================================================
T = table(all_group, all_infocap, all_suscep, ...
    'VariableNames', {'Group', 'Info_Cap', 'Susceptibility'});

% =========================================================
% Save CSV
% =========================================================
writetable(T, out_csv);

fprintf('CSV saved to:\n%s\n', out_csv);
fprintf('Total rows: %d\n', height(T));
disp(head(T));