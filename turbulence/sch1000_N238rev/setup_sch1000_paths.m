function P = setup_sch1000_paths()
%SETUP_SCH1000_PATHS Put active pipeline code ahead of historical copies.

P = pipeline_paths();

active_dirs = {
    P.sch1000_root
    P.prepare_data
    P.prepare_stage1
    P.prepare_stage2
    P.prepare_stage3
    P.prepare_stage4
    P.prepare_stage5
    P.calculate_turbu
    P.calculate_turbu_amyloid
    P.harmonization
    P.data_export
    P.neuromaps_export
    P.statistics
    P.visualization
    P.render
    P.render_surface
    P.render_utils
    P.helpers
    };

for i = 1:numel(active_dirs)
    if isfolder(active_dirs{i})
        addpath(active_dirs{i}, '-begin');
    end
end

if isfolder(P.render_utils)
    addpath(genpath(P.render_utils), '-begin');
end
if isfolder(P.render_assets) && ~strcmp(P.render_assets, P.render_utils)
    addpath(P.render_assets, '-begin');
end
% Rendering also contains historical copies of some renderer functions under
% render_utils. Put the maintained top-level implementations first.
if isfolder(P.render_surface)
    addpath(P.render_surface, '-begin');
end
end
