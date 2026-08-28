function demo_tbl = get_demographics_table(xlsx_file, sheet_name)
% LOAD_DEMOGRAPHICS_TABLE
%
% Expected columns:
%   PTID, GROUP, AGE, PTGENDER

    if nargin < 2 || isempty(sheet_name)
        demo_tbl = readtable(xlsx_file);
    else
        demo_tbl = readtable(xlsx_file, 'Sheet', sheet_name);
    end

    required_vars = {'PTID','GROUP','AGE','PTGENDER'};
    missing_vars = required_vars(~ismember(required_vars, demo_tbl.Properties.VariableNames));
    if ~isempty(missing_vars)
        error('Missing required variable(s): %s', strjoin(missing_vars, ', '));
    end

    demo_tbl = demo_tbl(:, required_vars);

    demo_tbl.PTID     = string(demo_tbl.PTID);
    demo_tbl.GROUP    = string(demo_tbl.GROUP);
    demo_tbl.PTGENDER = upper(strtrim(string(demo_tbl.PTGENDER)));

    demo_tbl = demo_tbl(~ismissing(demo_tbl.PTID) & ...
                        ~ismissing(demo_tbl.AGE) & ...
                        ~ismissing(demo_tbl.PTGENDER), :);

    % Numeric sex coding: Male = 1, Female = 0
    demo_tbl.SEX_NUM = nan(height(demo_tbl),1);
    demo_tbl.SEX_NUM(demo_tbl.PTGENDER == "M")      = 1;
    demo_tbl.SEX_NUM(demo_tbl.PTGENDER == "MALE")   = 1;
    demo_tbl.SEX_NUM(demo_tbl.PTGENDER == "F")      = 0;
    demo_tbl.SEX_NUM(demo_tbl.PTGENDER == "FEMALE") = 0;

    if any(isnan(demo_tbl.SEX_NUM))
        bad_vals = unique(demo_tbl.PTGENDER(isnan(demo_tbl.SEX_NUM)));
        error('Unrecognized PTGENDER value(s): %s', strjoin(cellstr(bad_vals), ', '));
    end
end