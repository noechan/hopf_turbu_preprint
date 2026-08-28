function paths = get_paths()

    paths.project_root = '/path/to/ADNI3/hopf_turbu_paper_project/main/';

    paths.code_root    = fullfile(paths.project_root, 'matlab/code');
    paths.data_root    = fullfile(paths.project_root, 'data', 'hopf_group');
    paths.precomp_root = fullfile(paths.data_root, 'precomputed_geometry');

    paths.results_root = fullfile(paths.project_root, 'results', 'hopf_group');
    paths.freq_root = fullfile(paths.project_root, 'results', 'hopf_group', 'freq');
    paths.empcorr_root = fullfile(paths.project_root, 'results', 'hopf_group', 'empirical_spacorr');
    paths.hopf_grange_root = fullfile(paths.project_root, 'results', 'hopf_group', 'G_range');
    paths.optG_root = fullfile(paths.project_root, 'results','hopf_group', 'optG');

    paths.abeta_table  = fullfile(paths.data_root, 'ADNI3_N238rev_with_ABETA_Status_CL24.xlsx');
    paths.precomp_file = fullfile(paths.precomp_root, 'precomputed_geometry.mat');

end