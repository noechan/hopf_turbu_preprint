function [age_all, sex_all, group_all] = get_covariates_for_ptids(ptid_list, demo_tbl)
% GET_COVARIATES_FOR_PTIDS
%
% Returns AGE, SEX_NUM and GROUP in the same order as ptid_list.

    ptid_list = string(ptid_list(:));
    demo_ptid = string(demo_tbl.PTID);

    [is_found, loc] = ismember(ptid_list, demo_ptid);

    if ~all(is_found)
        missing_ptids = unique(ptid_list(~is_found));
        error('PTIDs not found in demographics table: %s', ...
            strjoin(cellstr(missing_ptids), ', '));
    end

    age_all   = demo_tbl.AGE(loc);
    sex_all   = demo_tbl.SEX_NUM(loc);
    group_all = demo_tbl.GROUP(loc);
end