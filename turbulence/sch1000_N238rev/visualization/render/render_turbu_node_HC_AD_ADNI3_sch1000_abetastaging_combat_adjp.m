%% Manuscript Figure 4a: HC Aβ- versus AD Aβ+ node-wise map

clear all; close all; clc;

sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

pval_file = fullfile(P.nodewise_stats, ...
    'NodeWise_Turbu_lambda_0_01_HC_ABneg_vs_AD_ABpos_permFDR.xlsx');
output_path = fullfile(P.figures_harmonized, 'rendering_adjP');

assert(isfile(pval_file), 'Missing node-wise result: %s', pval_file);
if ~isfolder(output_path)
    mkdir(output_path);
end
cd(output_path);

S = readtable(pval_file);
assert(ismember('AdjP', S.Properties.VariableNames), ...
    'AdjP column missing from %s.', pval_file);

adj_p = S.AdjP(:);
assert(numel(adj_p) == 1000, 'Expected 1,000 adjusted p-values; found %d.', numel(adj_p));
assert(all(isfinite(adj_p) & adj_p > 0 & adj_p <= 1), ...
    'Adjusted p-values must be finite and in the interval (0, 1].');

% Figure 4 displays the complete statistical map, not a bottom-30% mask.
neglog_adjP = -log10(adj_p);
rangemin = min(neglog_adjP);
rangemax = max(neglog_adjP);

fprintf('HC- vs AD+, lambda=0.01: -log10(adjP) range %.6f to %.6f\n', ...
    rangemin, rangemax);

save('nodewise_neglog10_adjP_HC_ABneg_vs_AD_ABpos_lambda_0_01.mat', ...
    'adj_p', 'neglog_adjP', 'rangemin', 'rangemax', 'pval_file');

rendersurface_schaefer1000( ...
    neglog_adjP, rangemin, rangemax, 0, 'Purples9', 1);
print('render_HC_AD_adjP_neglog10_purples', '-dpng');
savefig('render_HC_AD_adjP_neglog10_purples.fig');
