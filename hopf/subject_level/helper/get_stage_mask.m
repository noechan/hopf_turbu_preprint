function [stage_mask, PTID] = get_stage_mask(cfg, paths)
% GET_STAGE_MASK
%
% Returns a logical mask selecting subjects belonging to the requested stage.

    abeta_table = readtable(paths.abeta_table);

    S = load(fullfile(paths.data_root, cfg.ptid_file));
    if ~isfield(S, 'PTID')
        error('Variable "PTID" not found in file: %s', fullfile(paths.data_root, cfg.ptid_file));
    end

    PTID = cellstr(S.PTID);

    stage_mask = ismember(PTID, ...
        abeta_table.PTID( ...
            strcmp(abeta_table.GROUP, cfg.group_label) & ...
            abeta_table.ABeta_pvc == cfg.abeta_value ...
        ) ...
    );
end