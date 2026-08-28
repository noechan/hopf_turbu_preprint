function P = pipeline_paths()
%PIPELINE_PATHS Resolve local code paths and mounted ADNI data paths.
%
% The repository may live anywhere. External ADNI data default to the
% mounted path below and can be overridden without editing code:
%
%   setenv('ADNI3_ROOT', '/path/to/ADNI3')

P.sch1000_root = fileparts(mfilename('fullpath'));
P.repo_root = fileparts(P.sch1000_root);

P.adni3_root = getenv('ADNI3_ROOT');
if isempty(P.adni3_root)
    P.adni3_root = fullfile(filesep, 'Volumes', 'ADNI', 'Projects', ...
        '2024', 'ADNI', 'LONI_IDA', 'ADNI3');
end

% Local, version-controlled code and manuscript-result locations.
P.prepare_data = fullfile(P.sch1000_root, 'prepare_data');
P.prepare_stage1 = fullfile(P.prepare_data, '01_cohort_selection');
P.prepare_stage2 = fullfile(P.prepare_data, '02_ptid_assembly');
P.prepare_stage3 = fullfile(P.prepare_data, '03_timeseries_extraction');
P.prepare_stage4 = fullfile(P.prepare_data, '04_combine_modalities');
P.prepare_stage5 = fullfile(P.prepare_data, '05_amyloid_groups');
P.calculate_turbu = fullfile(P.sch1000_root, 'calculate_turbu');
P.calculate_turbu_amyloid = fullfile(P.calculate_turbu, 'amyloid_status');
P.harmonization = fullfile(P.sch1000_root, 'harmonization_allfeat');
P.harmonization_results = fullfile(P.harmonization, 'recomputed');
P.data_export = fullfile(P.sch1000_root, 'data_export');
P.neuromaps_export = fullfile(P.data_export, 'neuromaps_exports');
P.statistics = fullfile(P.sch1000_root, 'statistical_analysis');
P.visualization = fullfile(P.sch1000_root, 'visualization');
P.render = fullfile(P.visualization, 'render');
P.render_surface = fullfile(P.render, 'RenderSurface');
P.figures = fullfile(P.sch1000_root, 'figures_N145');
P.figures_harmonized = fullfile(P.figures, 'sch1000', ...
    'Abeta_Status', 'harmonized_allfeat');
P.nodewise_stats = fullfile(P.figures_harmonized, 'nodewise_stats');
P.python_figures = fullfile(P.figures_harmonized, 'python');
P.figures_N238rev = fullfile(P.repo_root, 'figures_N238rev');
P.figures_matching = fullfile(P.repo_root, 'figures_HC_AD_MCI_all_matching');
P.helpers = fullfile(P.repo_root, 'helper_functions');

% External participant and time-series data on the ADNI drive.
P.participants = fullfile(P.adni3_root, 'participants');
P.participants_cross_sectional = fullfile(P.participants, 'cross_sectional');
P.preprocessed_timeseries = fullfile(P.adni3_root, 'timeseries', ...
    'prepro_fulldenoising');
P.timeseries_inputs = fullfile(P.adni3_root, 'timeseries', ...
    'inputs_fulldenoising');
P.timeseries_outputs = fullfile(P.adni3_root, 'timeseries', ...
    'outputs_fulldenoising');
P.harmonization_raw_inputs = fullfile(P.adni3_root, 'timeseries', ...
    'harmonization_inputs', 'sch1000_N238rev');
P.matching_sch1000_output = fullfile(P.timeseries_outputs, 'sch1000', ...
    'Turbu_by_node_ADNI3_matching_all_groups');
P.neuromaps_annotations = fullfile(P.adni3_root, 'neuromaps_analysis', ...
    'neuromaps-data', 'annotations', 'turbu', 'MNI152');
P.neuroharmonize_processed_turbu = fullfile(P.adni3_root, 'code', ...
    'ADNI3_neuroHarmonize_site', 'data', 'processed', 'turbu');
P.neuroharmonize_metadata = fullfile(P.adni3_root, 'code', ...
    'ADNI3_neuroHarmonize_site', 'data', 'raw', 'turbu');
P.prepare_metadata = fullfile(P.adni3_root, 'timeseries', ...
    'pipeline_metadata', 'sch1000_N238rev');

P.abeta_status_file = fullfile(P.participants_cross_sectional, ...
    'ADNI3_N238rev_with_ABETA_Status_CL24.xlsx');

P.input.HC = fullfile(P.timeseries_inputs, ...
    'ADNI3_HC_MPRAGE_IRFSPGR_N238rev');
P.input.MCI = fullfile(P.timeseries_inputs, ...
    'ADNI3_MCI_MPRAGE_IRFSPGR_N238rev');
P.input.AD = fullfile(P.timeseries_inputs, ...
    'ADNI3_AD_MPRAGE_IRFSPGR_N238rev');

P.output.HC = fullfile(P.timeseries_outputs, ...
    'Turbulence_ADNI3_HC_MPRAGE_IRFSPGR_N238rev_sch1000');
P.output.MCI = fullfile(P.timeseries_outputs, ...
    'Turbulence_ADNI3_MCI_MPRAGE_IRFSPGR_N238rev_sch1000');
P.output.AD = fullfile(P.timeseries_outputs, ...
    'Turbulence_ADNI3_AD_MPRAGE_IRFSPGR_N238rev_sch1000');

% Aβ-status analysis outputs. These are kept separate from the complete
% diagnosis-group outputs above because they contain different cohorts.
P.output.HC_ABetaNeg = fullfile(P.timeseries_outputs, ...
    'Turbulence_ADNI3_HC_ABetaNeg_MPRAGE_IRFSPGR_N238rev_sch1000');
P.output.HC_ABetaPos = fullfile(P.timeseries_outputs, ...
    'Turbulence_ADNI3_HC_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000');
P.output.MCI_ABetaPos = fullfile(P.timeseries_outputs, ...
    'Turbulence_ADNI3_MCI_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000');
P.output.AD_ABetaPos = fullfile(P.timeseries_outputs, ...
    'Turbulence_ADNI3_AD_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000');
end
