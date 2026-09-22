function hfig = rendersurface_dk68(dk68vector,rangemin,rangemax,inv,clmap,surfacetype,tick_format,render_assets)
%RENDERSURFACE_DK68 Render 68 cortical Desikan-Killiany parcel values.
%
% The vector order is 34 left-hemisphere parcels followed by the same 34
% right-hemisphere parcels. Within each hemisphere the order is the standard
% FreeSurfer aparc order used by abagen's Desikan-Killiany atlas; the corpus
% callosum label is excluded.

local_utils = fullfile(fileparts(mfilename('fullpath')), 'render_utils');
if ~exist('render_assets','var') || isempty(render_assets)
    render_assets = local_utils;
end
if ~exist('rangemin','var') || isempty(rangemin)
    rangemin = min(dk68vector);
end
if ~exist('rangemax','var') || isempty(rangemax)
    rangemax = max(dk68vector);
end
if ~exist('inv','var') || isempty(inv)
    inv = 0;
end
if ~exist('clmap','var') || isempty(clmap)
    clmap = 'Bu_10';
end
if ~exist('surfacetype','var') || isempty(surfacetype)
    surfacetype = 2;
end
if ~exist('tick_format','var')
    tick_format = '';
end

dk68vector = dk68vector(:);
assert(numel(dk68vector) == 68, ...
    'The Desikan-Killiany rendering vector must contain exactly 68 values.');
assert(all(isfinite(dk68vector)), 'All DK68 values must be finite.');
assert(isfinite(rangemin) && isfinite(rangemax) && rangemin < rangemax, ...
    'Color limits must be finite with rangemin < rangemax.');
assert(ismember(surfacetype, [1 2 3]), ...
    'surfacetype must be 1 (mid), 2 (inflated), or 3 (very inflated).');

required_assets = { ...
    'Glasser360.L.flat.32k_fs_LR.surf.gii', ...
    'Glasser360.L.inflated.32k_fs_LR.surf.gii', ...
    'Glasser360.L.mid.32k_fs_LR.surf.gii', ...
    'Glasser360.L.very_inflated.32k_fs_LR.surf.gii', ...
    'Glasser360.R.flat.32k_fs_LR.surf.gii', ...
    'Glasser360.R.inflated.32k_fs_LR.surf.gii', ...
    'Glasser360.R.mid.32k_fs_LR.surf.gii', ...
    'Glasser360.R.very_inflated.32k_fs_LR.surf.gii', ...
    'fsaverage.L.DKT_org_Atlas.32k_fs_LR.label.gii', ...
    'fsaverage.R.DKT_org_Atlas.32k_fs_LR.label.gii'};
missing_assets = required_assets(~cellfun( ...
    @(f) isfile(fullfile(render_assets, f)), required_assets));
if ~isempty(missing_assets)
    error(['Missing DK68 render assets in %s. Pass the asset directory as ' ...
        'the eighth input. First missing file: %s'], ...
        render_assets, missing_assets{1});
end
addpath(genpath(local_utils));
addpath(genpath(render_assets));

mid_l = gifti(fullfile(render_assets, ...
    'Glasser360.L.mid.32k_fs_LR.surf.gii'));
inflated_l = gifti(fullfile(render_assets, ...
    'Glasser360.L.inflated.32k_fs_LR.surf.gii'));
very_inflated_l = gifti(fullfile(render_assets, ...
    'Glasser360.L.very_inflated.32k_fs_LR.surf.gii'));
flat_l = gifti(fullfile(render_assets, ...
    'Glasser360.L.flat.32k_fs_LR.surf.gii'));
mid_r = gifti(fullfile(render_assets, ...
    'Glasser360.R.mid.32k_fs_LR.surf.gii'));
inflated_r = gifti(fullfile(render_assets, ...
    'Glasser360.R.inflated.32k_fs_LR.surf.gii'));
very_inflated_r = gifti(fullfile(render_assets, ...
    'Glasser360.R.very_inflated.32k_fs_LR.surf.gii'));
flat_r = gifti(fullfile(render_assets, ...
    'Glasser360.R.flat.32k_fs_LR.surf.gii'));

switch surfacetype
    case 1
        display_l = mid_l;
        display_r = mid_r;
    case 2
        display_l = inflated_l;
        display_r = inflated_r;
    case 3
        display_l = very_inflated_l;
        display_r = very_inflated_r;
end

label_l = gifti(fullfile(render_assets, ...
    'fsaverage.L.DKT_org_Atlas.32k_fs_LR.label.gii'));
label_r = gifti(fullfile(render_assets, ...
    'fsaverage.R.DKT_org_Atlas.32k_fs_LR.label.gii'));
values_l = nan(size(label_l.cdata));
values_r = nan(size(label_r.cdata));

% FreeSurfer labels 1:35 contain corpus callosum at 4, which is not one of
% the 34 cortical DK parcels returned by abagen.
dkt_labels = 1:35;
dkt_labels(4) = [];
for i = 1:34
    assert(any(label_l.cdata == dkt_labels(i)), ...
        'Left surface is missing DK label %d.', dkt_labels(i));
    assert(any(label_r.cdata == dkt_labels(i)), ...
        'Right surface is missing DK label %d.', dkt_labels(i));
    values_l(label_l.cdata == dkt_labels(i)) = dk68vector(i);
    values_r(label_r.cdata == dkt_labels(i)) = dk68vector(i + 34);
end

hfig = figure('Color', 'w', 'Position', [40, 40, 400, 600]);
layout = tiledlayout(hfig, 3, 2, 'TileSpacing', 'tight', ...
    'Padding', 'compact');
render_panel(nexttile(layout), display_l, values_l, [-90 0], ...
    rangemin, rangemax);
render_panel(nexttile(layout), display_l, values_l, [90 0], ...
    rangemin, rangemax);
render_panel(nexttile(layout), display_r, values_r, [-90 0], ...
    rangemin, rangemax);
render_panel(nexttile(layout), display_r, values_r, [90 0], ...
    rangemin, rangemax);
render_panel(nexttile(layout), flat_l, values_l, [0 90], ...
    rangemin, rangemax);
render_panel(nexttile(layout), flat_r, values_r, [0 90], ...
    rangemin, rangemax);

switch inv
    case 0
        if isnumeric(clmap)
            colors = clmap;
        else
            colors = othercolor(clmap);
        end
    case 1
        if isnumeric(clmap)
            colors = flipud(clmap);
        else
            colors = flipud(othercolor(clmap));
        end
    case 2
        assert(~isnumeric(clmap), ...
            'inv=2 requires a named othercolor colormap.');
        colors = othercolor(clmap, 3);
    otherwise
        error('inv must be 0, 1, or 2.');
end
colormap(hfig, colors);

cb = colorbar;
cb.Layout.Tile = 'south';
cb.Ticks = [rangemin rangemax];
if isempty(tick_format)
    cb.TickLabels = {sprintf('%.3f', rangemin), ...
        sprintf('%.3f', rangemax)};
else
    cb.TickLabels = {sprintf(tick_format, rangemin), ...
        sprintf(tick_format, rangemax)};
end
cb.FontSize = 14;
end


function render_panel(ax, surface, parcel_values, camera_view, rangemin, rangemax)
% Draw a neutral cortical surface, then overlay labeled parcel data.
axis(ax, 'equal');
axis(ax, 'off');
patch(ax, 'Faces', surface.faces, 'Vertices', surface.vertices, ...
    'FaceColor', [0.92 0.92 0.92], 'EdgeColor', 'none');
hold(ax, 'on');
patch(ax, 'Faces', surface.faces, 'Vertices', surface.vertices, ...
    'FaceVertexCData', parcel_values, 'FaceColor', 'interp', ...
    'EdgeColor', 'none');
hold(ax, 'off');
set(ax, 'CLim', [rangemin rangemax]);
view(ax, camera_view(1), camera_view(2));
camlight(ax);
lighting(ax, 'gouraud');
material(ax, 'dull');
end
