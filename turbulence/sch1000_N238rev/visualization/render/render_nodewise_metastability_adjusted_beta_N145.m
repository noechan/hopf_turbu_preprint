%% Covariate-adjusted node-level metastability differences (N = 145)
%
% This script renders the adjusted group coefficient from the current
% age-, sex-, and education-adjusted Freedman--Lane node-wise analyses.
% For a binary group indicator, Adjusted_Group_Beta is the covariate-
% adjusted mean difference Group_1 minus Group_2 at each Schaefer parcel.
%
% Both contrasts use the same colour limits so that their effect magnitudes
% can be compared directly. The complete 1,000-parcel vectors are rendered;
% no significance mask or top-30%% selection is applied.

clear all; close all; clc;

sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

assert(exist('gifti', 'file') ~= 0, ...
    ['The MATLAB GIFTI reader is not on the path. Add the SPM12 root or ' ...
     'another GIFTI toolbox before running this script.']);

results_dir = fullfile(P.statistics, 'nodewise_metastability', ...
    'results', 'N145_nodewise_metastability');
output_dir = fullfile(P.figures_harmonized, ...
    'rendering_adjusted_beta');

if ~isfolder(output_dir)
    mkdir(output_dir);
end

contrasts = struct( ...
    'id', { ...
        'HC_ABneg_vs_AD_ABpos', ...
        'MCI_ABpos_vs_AD_ABpos'}, ...
    'title', { ...
        'HC A\beta^- - AD A\beta^+', ...
        'MCI A\beta^+ - AD A\beta^+'}, ...
    'filename', { ...
        ['N145_nodewise_metastability_HC_ABneg_vs_AD_ABpos_' ...
         'FreedmanLane_age_sex_education.csv'], ...
        ['N145_nodewise_metastability_MCI_ABpos_vs_AD_ABpos_' ...
         'FreedmanLane_age_sex_education.csv']});

required_columns = { ...
    'Node', 'Group_1', 'Group_2', 'Adjusted_Group_Beta'};

for contrast_index = 1:numel(contrasts)
    input_file = fullfile(results_dir, contrasts(contrast_index).filename);
    assert(isfile(input_file), ...
        ['Missing node-wise result: %s\nRun ' ...
         'run_nodewise_metastability_N145_age_sex_education.R first.'], ...
        input_file);

    T = readtable(input_file);
    assert(all(ismember(required_columns, T.Properties.VariableNames)), ...
        'Missing required columns in %s.', input_file);
    assert(height(T) == 1000, ...
        'Expected 1,000 Schaefer parcels in %s; found %d.', ...
        input_file, height(T));

    node_names = string(T.Node);
    assert(all(startsWith(node_names, "Schaefer_")), ...
        'Unexpected parcel name in %s.', input_file);
    node_ids = str2double(extractAfter(node_names, "Schaefer_"));
    assert(all(isfinite(node_ids) & node_ids == round(node_ids)), ...
        'Parcel identifiers must be finite integers in %s.', input_file);
    assert(isequal(sort(node_ids), (1:1000)'), ...
        'Each Schaefer parcel from 1 through 1,000 must occur once in %s.', ...
        input_file);

    adjusted_beta = T.Adjusted_Group_Beta(:);
    assert(isnumeric(adjusted_beta) && all(isfinite(adjusted_beta)), ...
        'Adjusted_Group_Beta must contain 1,000 finite numeric values.');

    % Reorder explicitly by Schaefer parcel number. This prevents a changed
    % CSV row order from assigning values to the wrong surface parcels.
    beta_by_node = nan(1000, 1);
    beta_by_node(node_ids) = adjusted_beta;

    contrasts(contrast_index).input_file = input_file;
    contrasts(contrast_index).beta = beta_by_node;
    contrasts(contrast_index).group_1 = char(string(T.Group_1(1)));
    contrasts(contrast_index).group_2 = char(string(T.Group_2(1)));
end

all_beta = vertcat(contrasts.beta);

% Use one rounded scale for both contrasts. These manuscript contrasts are
% expected to be positive (earlier-stage group > AD+). If future reruns
% contain both signs, use a symmetric diverging scale automatically.
scale_step = 0.001;
if all(all_beta >= 0)
    rangemin = 0;
    rangemax = ceil(max(all_beta) / scale_step) * scale_step;
    colormap_name = 'Purples9';
    flip_colormap = 0;
else
    absolute_limit = ceil(max(abs(all_beta)) / scale_step) * scale_step;
    rangemin = -absolute_limit;
    rangemax = absolute_limit;
    colormap_name = 'RdBu11';
    % ColorBrewer RdBu is red-to-blue; flip it so positive effects are red.
    flip_colormap = 1;
end

surface_type = 2;       % inflated cortical surface
neutral_lowest = false; % zero/low effects are valid values, not a mask
tick_format = '%.3f';
render_assets = P.render_assets;

fprintf('Adjusted-beta common colour range: %.3f to %.3f\n', ...
    rangemin, rangemax);

for contrast_index = 1:numel(contrasts)
    adjusted_beta = contrasts(contrast_index).beta;
    fprintf('%s: adjusted-beta range %.6f to %.6f\n', ...
        contrasts(contrast_index).id, ...
        min(adjusted_beta), max(adjusted_beta));

    hfig = rendersurface_schaefer1000( ...
        adjusted_beta, rangemin, rangemax, flip_colormap, ...
        colormap_name, surface_type, neutral_lowest, tick_format, ...
        render_assets);
    set(hfig, 'Color', 'w');
    sgtitle(hfig, contrasts(contrast_index).title, ...
        'Interpreter', 'tex', 'FontWeight', 'bold', 'FontSize', 16);

    output_stem = fullfile(output_dir, sprintf( ...
        'N145_nodewise_metastability_adjusted_beta_%s', ...
        contrasts(contrast_index).id));
    print(hfig, [output_stem '.png'], '-dpng', '-r300');
    print(hfig, [output_stem '.pdf'], '-dpdf', '-r300', '-bestfit');
    savefig(hfig, [output_stem '.fig']);
end

metadata_file = fullfile(output_dir, ...
    'N145_nodewise_metastability_adjusted_beta_render_inputs.mat');
save(metadata_file, 'contrasts', 'rangemin', 'rangemax', ...
    'scale_step', 'surface_type', 'neutral_lowest', 'tick_format', ...
    'colormap_name', 'flip_colormap', 'render_assets', 'results_dir');

fprintf('Saved adjusted-difference renderings to: %s\n', output_dir);
fprintf('Saved rendering inputs and provenance to: %s\n', metadata_file);
