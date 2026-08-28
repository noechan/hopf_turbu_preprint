function report = preflight_sch1000()
%PREFLIGHT_SCH1000 Validate paths and core inputs without changing outputs.

P = setup_sch1000_paths();

checks = {
    'ADNI3 root', P.adni3_root, 'dir'
    'participant tables', P.participants_cross_sectional, 'dir'
    'preprocessed time series', P.preprocessed_timeseries, 'dir'
    'pipeline time-series inputs', P.timeseries_inputs, 'dir'
    'pipeline time-series outputs', P.timeseries_outputs, 'dir'
    'pipeline participant metadata', P.prepare_metadata, 'dir'
    'amyloid-status table', P.abeta_status_file, 'file'
    'Schaefer coordinates', fullfile(P.sch1000_root, 'SchaeferCOG.mat'), 'file'
    'Yeo network labels', fullfile(P.sch1000_root, 'RSN7vector.mat'), 'file'
    'HC Schaefer-1000 input', fullfile(P.input.HC, ...
        'tseries_ADNI3_HC_MPRAGE_IRFSPGR_sch1000_N238rev.mat'), 'file'
    'MCI Schaefer-1000 input', fullfile(P.input.MCI, ...
        'tseries_ADNI3_MCI_MPRAGE_IRFSPGR_sch1000_N238rev.mat'), 'file'
    'AD Schaefer-1000 input', fullfile(P.input.AD, ...
        'tseries_ADNI3_AD_MPRAGE_IRFSPGR_sch1000_N238rev.mat'), 'file'
    };

n = size(checks, 1);
ok = false(n, 1);
for i = 1:n
    if strcmp(checks{i, 3}, 'dir')
        ok(i) = isfolder(checks{i, 2});
    else
        ok(i) = isfile(checks{i, 2});
    end
    fprintf('[%s] %s\n       %s\n', pass_fail(ok(i)), ...
        checks{i, 1}, checks{i, 2});
end

required_functions = {'butter', 'filtfilt', 'hilbert', 'lsqcurvefit', ...
    'corr', 'ranksum'};
function_ok = false(numel(required_functions), 1);
for i = 1:numel(required_functions)
    function_ok(i) = ~isempty(which(required_functions{i}));
    fprintf('[%s] MATLAB function: %s\n', pass_fail(function_ok(i)), ...
        required_functions{i});
end

resolved_config = which('pipeline_paths');
local_config_ok = strcmp(resolved_config, fullfile(P.sch1000_root, ...
    'pipeline_paths.m'));
fprintf('[%s] Local cleaned repository is active\n       %s\n', ...
    pass_fail(local_config_ok), resolved_config);

report = struct();
report.Paths = table(string(checks(:, 1)), string(checks(:, 2)), ok, ...
    'VariableNames', {'Check', 'Path', 'Passed'});
report.Functions = table(string(required_functions(:)), function_ok, ...
    'VariableNames', {'Function', 'Passed'});
report.LocalCodeActive = local_config_ok;

if ~all(ok) || ~all(function_ok) || ~local_config_ok
    error('sch1000:PreflightFailed', ...
        'Preflight failed. Review the [FAIL] entries above.');
end

fprintf('\nSchaefer-1000 core preflight passed. No outputs were modified.\n');
end

function label = pass_fail(tf)
if tf
    label = 'PASS';
else
    label = 'FAIL';
end
end
