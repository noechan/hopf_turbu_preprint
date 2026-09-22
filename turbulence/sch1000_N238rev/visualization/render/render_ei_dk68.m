%% Transcriptomic E:I sanity-check maps in Desikan-Killiany 68
%
% Run build_ei_dk68.py before this script. The source map is calculated in
% the 34 left cortical DK parcels and reflected onto the homologous right
% parcels, matching the convention used in the reference paper.

clear all; close all; clc;

sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

assert(exist('gifti', 'file') ~= 0, ...
    ['The MATLAB GIFTI reader is not on the path. Add SPM or the GIFTI ' ...
     'toolbox before running this renderer.']);

ei_dir = fullfile(P.repository_root, 'abagen_analysis', 'abagen-code', ...
    'ei_maps', 'desikan_killiany_68');
ei_file = fullfile(ei_dir, 'EI_values_dk68_reflected.csv');
output_dir = fullfile(ei_dir, 'figures', 'surface_matlab');

assert(isfile(ei_file), ...
    'Missing DK68 E:I table. Run build_ei_dk68.py first: %s', ei_file);
if ~isfolder(output_dir)
    mkdir(output_dir);
end

T = readtable(ei_file);
required_columns = {'render_index', 'left_order', 'label', 'hemisphere', ...
    'source_left_label', 'ei_raw', 'ei_minmax'};
assert(all(ismember(required_columns, T.Properties.VariableNames)), ...
    'DK68 E:I table must contain columns: %s.', ...
    strjoin(required_columns, ', '));
assert(height(T) == 68, 'Expected 68 DK cortical parcels; found %d.', height(T));
assert(isequal(T.render_index(:), (1:68)'), ...
    'render_index must be ordered from 1 through 68.');
assert(all(string(T.hemisphere(1:34)) == "L") && ...
       all(string(T.hemisphere(35:68)) == "R"), ...
    'The first 34 rows must be left and the final 34 rows right hemisphere.');
assert(isequal(string(T.label(1:34)), string(T.label(35:68))), ...
    'Left and right DK parcel orders do not match.');

ei_raw = T.ei_raw(:);
ei_minmax = T.ei_minmax(:);
assert(all(isfinite(ei_raw)) && all(isfinite(ei_minmax)), ...
    'E:I values must all be finite.');
assert(max(abs(ei_raw(1:34) - ei_raw(35:68))) < 1e-12, ...
    'The right hemisphere is not an exact reflection of the left.');
assert(max(abs(ei_minmax(1:34) - ei_minmax(35:68))) < 1e-12, ...
    'The normalized right hemisphere is not an exact reflection of the left.');
assert(all(ei_minmax >= -eps & ei_minmax <= 1 + eps), ...
    'Min-max E:I values must lie within [0, 1].');

expected_minmax = (ei_raw - min(ei_raw)) ./ (max(ei_raw) - min(ei_raw));
assert(max(abs(ei_minmax - expected_minmax)) < 1e-10, ...
    'ei_minmax is not the expected min-max transform of ei_raw.');

raw_limits = [min(ei_raw), max(ei_raw)];
minmax_limits = [0, 1];
surface_type = 2;
flip_colormap = 0;
colormap_name = 'YlOrRd9';
render_assets = P.render_assets;

fprintf('Raw DK34 E:I range: %.6f to %.6f\n', ...
    raw_limits(1), raw_limits(2));

raw_figure = rendersurface_dk68(ei_raw, raw_limits(1), raw_limits(2), ...
    flip_colormap, colormap_name, surface_type, '%.3f', render_assets);
sgtitle(raw_figure, 'Raw transcriptomic E:I ratio (DK68 sanity check)', ...
    'FontWeight', 'bold', 'FontSize', 16);
raw_stem = fullfile(output_dir, 'EI_expression_raw_dk68_reflected_surface');
print(raw_figure, [raw_stem '.png'], '-dpng', '-r300');
print(raw_figure, [raw_stem '.pdf'], '-dpdf', '-r300', '-bestfit');
savefig(raw_figure, [raw_stem '.fig']);

minmax_figure = rendersurface_dk68(ei_minmax, 0, 1, ...
    flip_colormap, colormap_name, surface_type, '%.2f', render_assets);
sgtitle(minmax_figure, ...
    'Min-max transcriptomic E:I ratio (DK68 sanity check)', ...
    'FontWeight', 'bold', 'FontSize', 16);
minmax_stem = fullfile(output_dir, ...
    'EI_expression_minmax_dk68_reflected_surface');
print(minmax_figure, [minmax_stem '.png'], '-dpng', '-r300');
print(minmax_figure, [minmax_stem '.pdf'], '-dpdf', '-r300', '-bestfit');
savefig(minmax_figure, [minmax_stem '.fig']);

metadata_file = fullfile(output_dir, 'EI_surface_render_inputs.mat');
save(metadata_file, 'ei_raw', 'ei_minmax', 'raw_limits', ...
    'minmax_limits', 'surface_type', 'flip_colormap', 'colormap_name', ...
    'render_assets', 'ei_file');

fprintf('Saved DK68 E:I surface renderings to: %s\n', output_dir);
fprintf('Saved rendering inputs and provenance: %s\n', metadata_file);
