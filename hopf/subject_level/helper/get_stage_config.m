function cfg = get_stage_config(stage_key)
% GET_STAGE_CONFIG
%
% Returns configuration for one analysis stage.
%
% Valid stage_key values:
%   - 'HC_ABneg'
%   - 'HC_ABpos'
%   - 'MCI_ABpos'
%   - 'AD_ABpos'

    switch upper(stage_key)

        case 'HC_ABNEG'
            cfg.stage_key       = 'HC_ABneg';
            cfg.parent_group    = 'HC';
            cfg.group_label     = 'HC';
            cfg.abeta_value     = 0;   % Aβ−
            cfg.abeta_label     = 'ABneg';

            cfg.tseries_file    = 'tseries_ADNI3_HC_MPRAGE_IRFSPGR_sch1000_N238rev.mat';
            cfg.ptid_file       = 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat';
            cfg.gec_file    = 'ADNI3_HC_ABneg_GECglo_cnd_001_sch1000_mod_eps_1_2_ratio_SC_max02.mat';
            %cfg.lambda = 0.18;

            cfg.save_prefix     = 'ADNI3_HC_ABneg';
            cfg.results_subdir  = 'HC_ABneg';

        case 'HC_ABPOS'
            cfg.stage_key       = 'HC_ABpos';
            cfg.parent_group    = 'HC';
            cfg.group_label     = 'HC';
            cfg.abeta_value     = 1;   % Aβ+
            cfg.abeta_label     = 'ABpos';

            cfg.tseries_file    = 'tseries_ADNI3_HC_MPRAGE_IRFSPGR_sch1000_N238rev.mat';
            cfg.ptid_file       = 'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat';
            cfg.gec_file    = 'ADNI3_HC_ABpos_GECglo_cnd_001_sch1000_mod_eps_1_2_ratio_SC_max02.mat';
            %cfg.lambda = 0.18;

            cfg.save_prefix     = 'ADNI3_HC_ABpos';
            cfg.results_subdir  = 'HC_ABpos';

        case 'MCI_ABPOS'
            cfg.stage_key       = 'MCI_ABpos';
            cfg.parent_group    = 'MCI';
            cfg.group_label     = 'MCI';
            cfg.abeta_value     = 1;   % Aβ+
            cfg.abeta_label     = 'ABpos';

            cfg.tseries_file    = 'tseries_ADNI3_MCI_MPRAGE_IRFSPGR_sch1000_N238rev.mat';
            cfg.ptid_file       = 'PTID_ADNI3_MCI_MPRAGE_IRFSPGR_all.mat';
            cfg.gec_file    = 'ADNI3_MCI_ABpos_GECglo_cnd_001_sch1000_mod_eps_1_2_ratio_SC_max02.mat';

            cfg.save_prefix     = 'ADNI3_MCI_ABpos';
            cfg.results_subdir  = 'MCI_ABpos';

        case 'AD_ABPOS'
            cfg.stage_key       = 'AD_ABpos';
            cfg.parent_group    = 'AD';
            cfg.group_label     = 'AD';
            cfg.abeta_value     = 1;   % Aβ+
            cfg.abeta_label     = 'ABpos';

            cfg.tseries_file    = 'tseries_ADNI3_AD_MPRAGE_IRFSPGR_sch1000_N238rev.mat';
            cfg.ptid_file       = 'PTID_ADNI3_AD_MPRAGE_IRFSPGR_all.mat';
            cfg.gec_file    = 'ADNI3_AD_ABpos_GECglo_cnd_001_sch1000_mod_eps_1_2_ratio_SC_max02.mat';

            cfg.save_prefix     = 'ADNI3_AD_ABpos';
            cfg.results_subdir  = 'AD_ABpos';

        otherwise
            error('Unknown stage_key: %s', stage_key);
    end

    % Shared analysis parameters
    cfg.TR        = 3;
    cfg.NPARCELLS = 1000;
    cfg.Tmax      = 197;
    cfg.NCOND     = 1;
    cfg.NR        = 400;
    cfg.NRini     = 20;
    cfg.NRfin     = 380;
    cfg.NSUBSIM   = 100;
    cfg.lambda = 0.18;
end