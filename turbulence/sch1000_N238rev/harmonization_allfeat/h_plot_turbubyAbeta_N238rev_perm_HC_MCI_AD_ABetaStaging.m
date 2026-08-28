clear all; close all
rng(1, 'twister');  % Reproducible permutation p-values.
sch1000_root = fileparts(fileparts(mfilename('fullpath')));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();
output_path = fullfile(P.figures_harmonized, 'global');
if ~isfolder(output_path)
    mkdir(output_path);
end

%% Load results table
h_results_table = readtable(fullfile(P.harmonization_results, ...
    'Turbu_ComBat_ADNI3_allfeatures_N145.xlsx'));


%% ===== Extract variables from table for plotting =====
T = h_results_table;   % already loaded
assert(ismember('Group', T.Properties.VariableNames), 'Expected column "Group" not found.');

% --- Groups (exact labels provided) ---
idx_HC_ABneg  = strcmp(T.Group,'HC_ABneg');
idx_HC_ABpos  = strcmp(T.Group,'HC_ABpos');
idx_MCI_ABpos = strcmp(T.Group,'MCI_ABpos');
idx_AD_ABpos  = strcmp(T.Group,'AD_ABpos');

% --- Lambda definitions (titles use dot; headers use underscore) ---
lambda_dot = {'0.27','0.24','0.21','0.18','0.15','0.12','0.09','0.06','0.03','0.01'};
lambda_us  = strrep(lambda_dot,'.','_');  % {'0_27', '0_24', ...}

% --- Helper to pull a lambda-resolved matrix (rows = lambdas, cols = subjects) ---
getMat = @(prefix, mask, lam_us) ...
    ( T{mask, strcat(prefix, lam_us)}.' );

% ----- Turbulence -----
Turbu_HC_ABneg   = getMat('Turbu_lam_',       idx_HC_ABneg,  lambda_us);
Turbu_HC_ABpos   = getMat('Turbu_lam_',       idx_HC_ABpos,  lambda_us);
Turbu_MCI_ABpos  = getMat('Turbu_lam_',       idx_MCI_ABpos, lambda_us);
Turbu_AD_ABpos   = getMat('Turbu_lam_',       idx_AD_ABpos,  lambda_us);

% ----- Information Flow -----
lam_if_idx = 2:numel(lambda_us);                    % skip 0.27
Info_flow_HC_ABneg   = getMat('InfoFlow_lam_',     idx_HC_ABneg,  lambda_us(lam_if_idx));
Info_flow_HC_ABpos   = getMat('InfoFlow_lam_',     idx_HC_ABpos,  lambda_us(lam_if_idx));
Info_flow_MCI_ABpos  = getMat('InfoFlow_lam_',     idx_MCI_ABpos, lambda_us(lam_if_idx));
Info_flow_AD_ABpos   = getMat('InfoFlow_lam_',     idx_AD_ABpos,  lambda_us(lam_if_idx));
% Back-compat alias used later in your script:
Info_flow_AD = Info_flow_AD_ABpos;

% ----- Information Transfer (includes 0.27 in your headers) -----
Info_transfer_HC_ABneg   = getMat('InfoTransfer_lam_', idx_HC_ABneg,  lambda_us);
Info_transfer_HC_ABpos   = getMat('InfoTransfer_lam_', idx_HC_ABpos,  lambda_us);
Info_transfer_MCI_ABpos  = getMat('InfoTransfer_lam_', idx_MCI_ABpos, lambda_us);
Info_transfer_AD_ABpos   = getMat('InfoTransfer_lam_', idx_AD_ABpos,  lambda_us);
% Back-compat alias:
Info_transfer_AD = Info_transfer_AD_ABpos;

% ----- Scalars -----
Info_cascade_HC_ABneg  = T.InfoCascade(idx_HC_ABneg);
Info_cascade_HC_ABpos  = T.InfoCascade(idx_HC_ABpos);
Info_cascade_MCI_ABpos = T.InfoCascade(idx_MCI_ABpos);
Info_cascade_AD_ABpos  = T.InfoCascade(idx_AD_ABpos);

Metastability_HC_ABneg  = T.Metastability(idx_HC_ABneg);
Metastability_HC_ABpos  = T.Metastability(idx_HC_ABpos);
Metastability_MCI_ABpos = T.Metastability(idx_MCI_ABpos);
Metastability_AD_ABpos  = T.Metastability(idx_AD_ABpos); 

gKoP_HC_ABneg  = T.gKoP(idx_HC_ABneg);
gKoP_HC_ABpos  = T.gKoP(idx_HC_ABpos);
gKoP_MCI_ABpos = T.gKoP(idx_MCI_ABpos);
gKoP_AD_ABpos  = T.gKoP(idx_AD_ABpos);                 

% --- Make sure the "lambda" cell array in your plotting section matches the titles ---
lambda = lambda_dot;   % reuse for subplot titles and legends

% (Optional) quick sanity printout
fprintf('N per group — HC-: %d | HC+: %d | MCI+: %d | AD+: %d\n', ...
    nnz(idx_HC_ABneg), nnz(idx_HC_ABpos), nnz(idx_MCI_ABpos), nnz(idx_AD_ABpos));

cd(output_path)

%% Plotting setup
lambda = {'0.27','0.24','0.21','0.18','0.15','0.12','0.09','0.06','0.03','0.01'};
subgroups4 = {'HC-AB-','HC-AB+','MCI-AB+','AD-AB+'};

% Base-MATLAB padding; avoids the Image Processing Toolbox dependency of padarray.
padToMat = @(C) cell2mat(cellfun(@(x) ...
    [x(:); NaN(max(cellfun(@numel, C)) - numel(x), 1)], ...
    C, 'UniformOutput', false));

%% --- Turbulence (4 groups) ---
figure
p1P = nan(5, numel(lambda)); 
for lamb = 1:numel(lambda)
    subplot(2,5,lamb)

    C = {Turbu_HC_ABneg(lamb,:), Turbu_HC_ABpos(lamb,:), Turbu_MCI_ABpos(lamb,:), Turbu_AD_ABpos(lamb,:)};
    Cmat = padToMat(C);

    boxplot(Cmat, 'Labels', subgroups4, 'Symbol',''); title(lambda{lamb});

    [p_HCneg_HCpos,~,~] = permutationTest(Turbu_HC_ABneg(lamb,:), Turbu_HC_ABpos(lamb,:), 1000);
    [p_HCpos_MCIpos,~,~]= permutationTest(Turbu_HC_ABpos(lamb,:), Turbu_MCI_ABpos(lamb,:), 1000);
    [p_MCIpos_ADpos,~,~] = permutationTest(Turbu_MCI_ABpos(lamb,:), Turbu_AD_ABpos(lamb,:), 1000);
    [p_HCneg_MCIpos,~,~] = permutationTest(Turbu_HC_ABneg(lamb,:), Turbu_MCI_ABpos(lamb,:), 1000);
    [p_HCneg_ADpos,~,~] = permutationTest(Turbu_HC_ABneg(lamb,:), Turbu_AD_ABpos(lamb,:), 1000);

    p1P(:,lamb) = [p_HCneg_HCpos; p_HCpos_MCIpos; p_MCIpos_ADpos;p_HCneg_MCIpos;p_HCneg_ADpos];

    sigstar({[1,2],[2,3],[3,4],[1,3],[1,4]}, ...
            [p_HCneg_HCpos,p_HCpos_MCIpos,p_MCIpos_ADpos,p_HCneg_MCIpos,p_HCneg_ADpos]);
end
print('Turbulence_HC_MCI_AD_ABetaStaging_Permutation 5 tests','-dpng')
saveas(gcf,'Turbulence_HC_MCI_AD_ABetaStaging_Permutation 5 tests','fig')
close

% Example FDR across selected lambdas (adjust as needed)
sel = [1:10]; % 0.27, 0.24, 0.09, 0.01
[h_t,~,~,adj_p_t] = fdr_bh(p1P(1,sel),0.05,'pdep','yes'); % HC- vs HC+
save('Turbu_Global_rejectedH0_HCabneg_vs_HCabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p1P(2,sel),0.05,'pdep','yes'); % HC+ vs MCI+
save('Turbu_Global_rejectedH0_HCabpos_vs_MCIabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p1P(3,sel),0.05,'pdep','yes'); % MCI+ vs AD+
save('Turbu_Global_rejectedH0_MCIabpos_vs_ADabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p1P(4,sel),0.05,'pdep','yes'); % HC- vs MCI+
save('Turbu_Global_rejectedH0_HCabneg_vs_MCIabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p1P(5,sel),0.05,'pdep','yes'); % HC- vs AD+
save('Turbu_Global_rejectedH0_HCabneg_vs_ADabpos.mat','h_t','adj_p_t')

%% --- Information Flow (lambda-resolved) ---
figure
p2P = nan(5, numel(lambda)-1);
for lamb = 2:numel(lambda)
    row = lamb - 1;     
    col = lamb - 1;    
    subplot(2,5,lamb)

    C = {Info_flow_HC_ABneg(row,:), Info_flow_HC_ABpos(row,:), Info_flow_MCI_ABpos(row,:), Info_flow_AD_ABpos(row,:)};
    Cmat = padToMat(C);

    boxplot(Cmat, 'Labels', subgroups4, 'Symbol',''); title(lambda{lamb});

    [p_HCneg_HCpos,~,~] = permutationTest(Info_flow_HC_ABneg(row,:), Info_flow_HC_ABpos(row,:), 10000);
    [p_HCpos_MCIpos,~,~]= permutationTest(Info_flow_HC_ABpos(row,:), Info_flow_MCI_ABpos(row,:), 10000);
    [p_MCIpos_AD,~,~]   = permutationTest(Info_flow_MCI_ABpos(row,:), Info_flow_AD(row,:), 10000);
    [p_HCneg_MCIpos,~,~] = permutationTest(Info_flow_HC_ABneg(row,:), Info_flow_MCI_ABpos(row,:), 10000);
    [p_HCneg_ADpos,~,~] = permutationTest(Info_flow_HC_ABneg(row,:), Info_flow_AD_ABpos(row,:), 10000);

    p2P(:,col) = [p_HCneg_HCpos;p_HCpos_MCIpos; p_MCIpos_AD;p_HCneg_MCIpos;p_HCneg_ADpos];

    sigstar({[1,2],[2,3],[3,4],[1,3],[1,4]}, ...
            [p_HCneg_HCpos,p_HCpos_MCIpos,p_MCIpos_AD,p_HCneg_MCIpos,p_HCneg_ADpos]);
end
print('InfoFlow_HC_MCI_AD_ABetaStaging_Permutation 5 tests','-dpng')
saveas(gcf,'InfoFlow_HC_MCI_AD_ABetaStaging_Permutation 5 tests','fig')
close

%FDR correction
sel = 1:9;  
[h_t,~,~,adj_p_t] = fdr_bh(p2P(1,sel),0.05,'pdep','yes'); % HC- vs HC+
save('Info_flow_rejectedH0_HCabneg_vs_HCabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p2P(2,sel),0.05,'pdep','yes'); % HC+ vs MCI+
save('Info_flow_rejectedH0_HCabpos_vs_MCIabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p2P(3,sel),0.05,'pdep','yes'); % MCI+ vs AD+
save('Info_flow_rejectedH0_MCIabpos_vs_ADabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p2P(4,sel),0.05,'pdep','yes'); % HC- vs MCI+
save('Info_flow_rejectedH0_HCabneg_vs_MCIabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p2P(5,sel),0.05,'pdep','yes'); % HC- vs AD+
save('Info_flow_rejectedH0_HCabneg_vs_ADabpos.mat','h_t','adj_p_t')

% Shaded mean±SD curves
LAMBDA = [0.24,0.21,0.18,0.15,0.12,0.09,0.06,0.03,0.01];
figure;
shadedErrorBar(LAMBDA(:), Info_flow_HC_ABneg(:,:)', {@mean,@std}, 'lineprops',{'-o','markerfacecolor','auto'}); hold on
shadedErrorBar(LAMBDA(:), Info_flow_HC_ABpos(:,:)', {@mean,@std}, 'lineprops',{'-o','markerfacecolor','auto'});
shadedErrorBar(LAMBDA(:), Info_flow_MCI_ABpos(:,:)', {@mean,@std}, 'lineprops',{'-o','markerfacecolor','auto'});
shadedErrorBar(LAMBDA(:), Info_flow_AD_ABpos(:,:)',       {@mean,@std}, 'lineprops',{'-o','markerfacecolor','auto'});
ylabel('Info Flow'); xlabel('Lambda Values')
legend({'HC-AB-','HC-AB+','MCI-AB+','AD'}, 'location','best')
print('InfoFlow_HC_MCI_AD_ABetaStaging_Permutation_ShadedErrorBar','-dpng')
saveas(gcf,'InfoFlow_HC_MCI_AD_ABetaStaging_Permutation_ShadedErrorBar','fig')
close

%% --- Information Cascade (scalar) ---
figure
C = {Info_cascade_HC_ABneg, Info_cascade_HC_ABpos, Info_cascade_MCI_ABpos, Info_cascade_AD_ABpos};
Cmat = padToMat(C);
boxplot(Cmat, 'Labels', subgroups4, 'Symbol','');
title('Information Cascade by ABETA Status')

[p_HCneg_HCpos,~,~] = permutationTest(Info_cascade_HC_ABneg, Info_cascade_HC_ABpos, 10000);
[p_HCpos_MCIpos,~,~]= permutationTest(Info_cascade_HC_ABpos, Info_cascade_MCI_ABpos, 10000);
[p_MCIpos_AD,~,~]   = permutationTest(Info_cascade_MCI_ABpos, Info_cascade_AD_ABpos, 10000);
[p_HCneg_MCIpos,~,~] = permutationTest(Info_cascade_HC_ABneg, Info_cascade_MCI_ABpos, 10000);
[p_HCneg_ADpos,~,~] = permutationTest(Info_cascade_HC_ABneg, Info_cascade_AD_ABpos, 10000);

p3P = [p_HCneg_HCpos; p_HCpos_MCIpos; p_MCIpos_AD;p_HCneg_MCIpos;p_HCneg_ADpos];
sigstar({[1,2],[2,3],[3,4],[1,3],[1,4]}, p3P');

print('Info_Cascade_HC_MCI_AD_ABetaStaging_Permutation 5 tests','-dpng')
saveas(gcf,'Info_Cascade_HC_MCI_AD_ABetaStaging_Permutation 5 tests','fig')
close

%% --- Information Transfer (lambda-resolved; plot 1 - value as in your code) ---
figure
p4P = nan(5, numel(lambda));
for lamb = 1:numel(lambda)
    subplot(2,5,lamb)

    C = {1-Info_transfer_HC_ABneg(lamb,:), 1-Info_transfer_HC_ABpos(lamb,:), 1-Info_transfer_MCI_ABpos(lamb,:), 1-Info_transfer_AD_ABpos(lamb,:)};
    Cmat = padToMat(C);

    boxplot(Cmat, 'Labels', subgroups4, 'Symbol',''); title(lambda{lamb});

    [p_HCneg_HCpos,~,~] = permutationTest(Info_transfer_HC_ABneg(lamb,:), Info_transfer_HC_ABpos(lamb,:), 1000);
    [p_HCpos_MCIpos,~,~]= permutationTest(Info_transfer_HC_ABpos(lamb,:), Info_transfer_MCI_ABpos(lamb,:), 1000);
    [p_MCIpos_AD,~,~]   = permutationTest(Info_transfer_MCI_ABpos(lamb,:), Info_transfer_AD(lamb,:), 1000);
    [p_HCneg_MCIpos,~,~] = permutationTest(Info_transfer_HC_ABneg(lamb,:), Info_transfer_MCI_ABpos(lamb,:), 10000);
    [p_HCneg_ADpos,~,~] = permutationTest(Info_transfer_HC_ABneg(lamb,:), Info_transfer_AD_ABpos(lamb,:), 10000);
    p4P(:,lamb) = [p_HCneg_HCpos; p_HCpos_MCIpos; p_MCIpos_AD;p_HCneg_MCIpos;p_HCneg_ADpos];

    sigstar({[1,2],[2,3],[3,4],[1,3],[1,4]}, ...
            [p_HCneg_HCpos,p_HCpos_MCIpos,p_MCIpos_AD,p_HCneg_MCIpos,p_HCneg_ADpos]);
end
print('InfoTransfer_HC_MCI_AD_ABetaStaging_Permutation 5 tests','-dpng')
saveas(gcf,'InfoTransfer_HC_MCI_AD_ABetaStaging_Permutation 5 tests','fig')
close

%FDR correction
sel = [1:10]; 
[h_t,~,~,adj_p_t] = fdr_bh(p4P(1,sel),0.05,'pdep','yes'); % HC- vs HC+
save('Info_transfer_rejectedH0_HCabneg_vs_HCabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p4P(2,sel),0.05,'pdep','yes'); % HC+ vs MCI+
save('Info_transfer_rejectedH0_HCabpos_vs_MCIabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p4P(3,sel),0.05,'pdep','yes'); % MCI+ vs AD+
save('Info_transfer_rejectedH0_MCIabpos_vs_ADabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p4P(4,sel),0.05,'pdep','yes'); % HC- vs MCI+
save('Info_transfer_rejectedH0_HCabneg_vs_MCIabpos.mat','h_t','adj_p_t')
[h_t,~,~,adj_p_t] = fdr_bh(p4P(5,sel),0.05,'pdep','yes'); % HC- vs AD+
save('Info_transfer_rejectedH0_HCabneg_vs_ADabpos.mat','h_t','adj_p_t')

%% --- Metastability (scalar) ---
figure
C = {Metastability_HC_ABneg, Metastability_HC_ABpos, Metastability_MCI_ABpos, Metastability_AD_ABpos};
Cmat = padToMat(C);
boxplot(Cmat, 'Labels', subgroups4, 'Symbol',''); title('Metastability by ABETA Status')

[p_HCneg_HCpos,~,~] = permutationTest(Metastability_HC_ABneg, Metastability_HC_ABpos, 1000);
[p_HCpos_MCIpos,~,~]= permutationTest(Metastability_HC_ABpos, Metastability_MCI_ABpos, 1000);
[p_MCIpos_ADpos,~,~]   = permutationTest(Metastability_MCI_ABpos, Metastability_AD_ABpos, 1000);

p5P = [p_HCneg_HCpos; p_HCpos_MCIpos; p_MCIpos_ADpos];
sigstar({[1,2],[2,3],[3,4]}, p5P');

print('Metastability_HCpm_AD_MCIabpos_Permutation 5 tests','-dpng')
saveas(gcf,'Metastability_HCpm_AD_MCIabpos_Permutation 5 tests','fig')
close

%% --- Global Kuramoto (gKoP; scalar) ---
figure
C = {gKoP_HC_ABneg, gKoP_HC_ABpos, gKoP_MCI_ABpos, gKoP_AD_ABpos};
Cmat = padToMat(C);
boxplot(Cmat, 'Labels', subgroups4, 'Symbol',''); title('Global Kuramoto by ABETA Status')

[p_HCneg_HCpos,~,~] = permutationTest(gKoP_HC_ABneg, gKoP_HC_ABpos, 1000);
[p_HCpos_MCIpos,~,~]= permutationTest(gKoP_HC_ABpos, gKoP_MCI_ABpos, 1000);
[p_MCIpos_AD,~,~]   = permutationTest(gKoP_MCI_ABpos, gKoP_AD_ABpos, 1000);

p6P = [p_HCneg_HCpos;p_HCpos_MCIpos; p_MCIpos_AD];
sigstar({[1,2],[2,3],[3,4]}, p6P');

print('gKoP_HCpm_AD_MCIabpos_Permutation 5 tests','-dpng')
saveas(gcf,'gKoP_HCpm_AD_MCIabpos_Permutation 5 tests','fig')
close
