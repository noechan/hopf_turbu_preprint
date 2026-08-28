%% ==========================================================
%  Build PTID + CAP + Susceptibility tables (separate folders)
% ==========================================================

clear; clc;

outDir = fullfile(pwd,'ML');
if ~exist(outDir,'dir'); mkdir(outDir); end

% ---- Define file locations ----
groups = struct([]);

groups(1).name = "HC_neg";
groups(1).pert_file = "/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/plots/perturbation_HC_ABneg_SUB.mat";
groups(1).emp_file  = "/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/data/empirical_spacorr_rest_cond_1_ADNI3_HC_ABneg_sch1000_SUB.mat";

groups(2).name = "HC_pos";
groups(2).pert_file = "/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/plots/perturbation_HC_ABpos_SUB.mat";
groups(2).emp_file  = "/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/data/empirical_spacorr_rest_cond_1_ADNI3_HC_ABpos_sch1000_SUB.mat";

groups(3).name = "MCI_pos";
groups(3).pert_file = "/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/plots/perturbation_MCI_ABpos_SUB.mat";
groups(3).emp_file  = "/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/data/empirical_spacorr_rest_cond_1_ADNI3_MCI_ABpos_sch1000_SUB.mat";

groups(4).name = "AD_pos";
groups(4).pert_file = "/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/plots/perturbation_AD_ABpos_SUB.mat";
groups(4).emp_file  = "/path/to/ADNI3/code/HPC_Hopf_SUB_DTI_1000_Staging/data/empirical_spacorr_rest_cond_1_ADNI3_AD_ABpos_sch1000_SUB.mat";


% ---- Run ----
for g = 1:numel(groups)

    fprintf('\nProcessing %s...\n',groups(g).name)

    T = build_cap_sus_table(groups(g).pert_file, groups(g).emp_file);

    save(fullfile(outDir, groups(g).name+"_CAP_Sus.mat"),'T','-v7.3');
    excelFile = fullfile(outDir, groups(g).name + "_CAP_Sus.xlsx");
    writetable(T, excelFile, 'FileType','spreadsheet', 'Sheet','CAP_Sus');


end

fprintf('\nAll groups finished.\n');
