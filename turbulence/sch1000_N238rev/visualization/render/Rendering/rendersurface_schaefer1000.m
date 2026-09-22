function hfig = rendersurface_schaefer1000(schaefer1000vector,rangemin,rangemax,inv,clmap,surfacetype,neutral_lowest,tick_format,render_assets)
% script for rendering a schaefer1000vector
%   ML Kringelbach April 2020
% schaefer1000vector : the schaefer1000vector values to be rendered
%
% rangemin, rangemax : limits for colorscheme
%
% inv:
%  0 colormap, interp
%  1 flip colormap, interp
%  2 colormap, only three colours
%
% clmap : name accepted by othercolor (default 'Bu_10') or an N-by-3 RGB
% matrix such as parula(256)
%
% surfacetype: 
%  1 midthickness
%  2 inflated (default)
%  3 very inflated
%
% neutral_lowest:
%  true replaces the lowest colormap entry with light grey (default,
%  retained for compatibility with the turbulence statistical maps)
%  false preserves the complete sequential colormap
%
% tick_format:
%  optional sprintf format applied to both colorbar endpoints. An empty
%  value retains the historical %.0f minimum / %.3f maximum formatting.
%
% render_assets:
%  directory containing the surface and Schaefer label GIFTI files.
%  Defaults to RenderSurface/render_utils for backward compatibility.

% The publication repository bundles the surface assets under render_utils.
local_utils = fullfile(fileparts(mfilename('fullpath')), 'render_utils');
if ~exist('render_assets','var') || isempty(render_assets)
    render_assets = local_utils;
end
render_utils = render_assets;
required_assets = { ...
    'Glasser360.L.flat.32k_fs_LR.surf.gii', ...
    'Glasser360.L.inflated.32k_fs_LR.surf.gii', ...
    'Glasser360.L.mid.32k_fs_LR.surf.gii', ...
    'Glasser360.L.very_inflated.32k_fs_LR.surf.gii', ...
    'Glasser360.R.flat.32k_fs_LR.surf.gii', ...
    'Glasser360.R.inflated.32k_fs_LR.surf.gii', ...
    'Glasser360.R.mid.32k_fs_LR.surf.gii', ...
    'Glasser360.R.very_inflated.32k_fs_LR.surf.gii', ...
    'Schaefer1000_L.func.gii', 'Schaefer1000_R.func.gii'};
missing_assets = required_assets(~cellfun(@(f) isfile(fullfile(render_utils, f)), required_assets));
if ~isempty(missing_assets)
    error(['Missing Schaefer render assets in %s. Pass the asset directory ' ...
        'as the ninth input. First missing file: %s'], ...
        render_utils, missing_assets{1});
end
addpath(genpath(local_utils))
addpath(genpath(render_utils))


if ~exist('rangemin','var')
     rangemin=min(schaefer1000vector);
end

if ~exist('rangemax','var')
     rangemax=max(schaefer1000vector);
end

if ~exist('inv','var')
     inv=max(schaefer1000vector);
end

if ~exist('clmap','var')
    clmap='Bu_10';
end

if ~exist('surfacetype','var')
     surfacetype=2; % default is inflated
end

if ~exist('neutral_lowest','var') || isempty(neutral_lowest)
    neutral_lowest=true;
end

if ~exist('tick_format','var')
    tick_format='';
end


% make space tight
make_it_tight = true;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
if ~make_it_tight,  clear subplot;  end

basedir = render_utils;
glassers_L=gifti(fullfile(basedir, 'Glasser360.L.mid.32k_fs_LR.surf.gii'));
glassersi_L=gifti(fullfile(basedir, 'Glasser360.L.inflated.32k_fs_LR.surf.gii'));
glassersvi_L=gifti(fullfile(basedir, 'Glasser360.L.very_inflated.32k_fs_LR.surf.gii'));
glassersf_L=gifti(fullfile(basedir, 'Glasser360.L.flat.32k_fs_LR.surf.gii'));
glassers_R=gifti(fullfile(basedir, 'Glasser360.R.mid.32k_fs_LR.surf.gii'));
glassersi_R=gifti(fullfile(basedir, 'Glasser360.R.inflated.32k_fs_LR.surf.gii'));
glassersvi_R=gifti(fullfile(basedir, 'Glasser360.R.very_inflated.32k_fs_LR.surf.gii'));
glassersf_R=gifti(fullfile(basedir, 'Glasser360.R.flat.32k_fs_LR.surf.gii'));

switch surfacetype
    case 1
        display_surf_left=glassers_L;
        display_surf_right=glassers_R;
    case 2
        display_surf_left=glassersi_L;
        display_surf_right=glassersi_R;
    case 3
        display_surf_left=glassersvi_L;
        display_surf_right=glassersvi_R;
end

sl = display_surf_left;
sr = display_surf_right;


% base='/path/to/local/schaefer/resources/';
base = render_utils;
atlas_l=gifti(fullfile(base, 'Schaefer1000_L.func.gii'));
atlas_r=gifti(fullfile(base, 'Schaefer1000_R.func.gii'));

% The former renderer loaded a Schaefer CIFTI label file here using
% FieldTrip's ft_read_cifti, but never used the returned value. Parcel
% labels are already supplied by the left/right GIFTI atlas files above.
% Omitting that dead read removes an unnecessary FieldTrip dependency.

vl=atlas_l;
vr=atlas_r;


% % replace the left hemisphere labels with values from dkt(1:34)
% % dk46_L.labels contains 34 elements
% left hemisphere
for i=1:500
    idx=find(atlas_l.cdata==i);
    vl.cdata(idx)=schaefer1000vector(i);
end

% right hemisphere
for i=501:1000
    idx=find(atlas_r.cdata==i);
    vr.cdata(idx)=schaefer1000vector(i);
end

% remove -1 from labels
% idx{i}=find(label_L.cdata==-1);
% vl.cdata(idx{i})=0;
% idx{i}=find(label_R.cdata==-1);
% vr.cdata(idx{i})=0;


%% rendering

    % create figure
    hfig = figure;
    set(gcf, 'Position',  [40, 40, 400, 600]);
    t=tiledlayout(3,2);
    t.TileSpacing='tight';

    %subplot(3,2,1); %left hemisphere side view
    nexttile
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(1) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(-90,0);
    camlight;
    lighting gouraud;
    material dull;


    %subplot(3,2,3); %left hemisphere midline
    nexttile
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(1) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(90,0)
    camlight;
    lighting gouraud;
    material dull;

    %subplot(3,2,4); %right hemisphere side view
    nexttile
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(-90,0)
    camlight;
    lighting gouraud;
    material dull;
 
    %subplot(3,2,2); %right hemisphere midline
    nexttile
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(90,0)
    camlight;
    lighting gouraud;
    material dull;

    % flatmaps
    sl=glassersf_L;
    sr=glassersf_R;
    
    %subplot(3,2,5); %left hemisphere flat
    nexttile
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sl.faces,'vertices',sl.vertices, 'FaceVertexCData', vl.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(0,90)
    camlight;
    lighting gouraud;
    material dull;
    %colorbar('southoutside');

    %subplot(3,2,6); %right hemisphere flat
    nexttile
    ax2=gca;
    axis(ax2,'equal');
    axis(ax2,'off');
    s(2) = patch(ax2,'Faces',sr.faces,'vertices',sr.vertices, 'FaceVertexCData', vr.cdata, 'FaceColor','interp', 'EdgeColor', 'none');
    set(ax2,'CLim',[rangemin rangemax]);
    view(0,90)
    camlight;
    lighting gouraud;
    material dull;
    %colorbar('southoutside');

    cb=colorbar;
    cb.Layout.Tile = 3;
    cb.Layout.Tile = 'south';
    cb.Ticks=([rangemin rangemax]);
    if isempty(tick_format)
        tmin=sprintf('%.0f',rangemin);
        tmax=sprintf('%.3f',rangemax);
    else
        tmin=sprintf(tick_format,rangemin);
        tmax=sprintf(tick_format,rangemax);
    end
    cb.TickLabels={tmin,tmax};
    cb.FontSize=14;

    switch inv
        case 0
            % use selected colormap (interpolated to 64 values)
            if isnumeric(clmap)
                c=clmap;
            else
                c=othercolor(clmap);
            end
        case 1
            % flip colormap (interpolated to 64 values)
            if isnumeric(clmap)
                c=flipud(clmap);
            else
                c=flipud(othercolor(clmap));
            end
        case 2
            % use specialised version with only three values
            assert(~isnumeric(clmap), ...
                'inv=2 requires a named othercolor colormap.');
            c=othercolor(clmap,3);
            % this is for the second value
            c(2,1)=0.70;c(2,2)=0.70;c(2,3)=0.70; 
        case 3
            % use specialised version with only four values
            assert(~isnumeric(clmap), ...
                'inv=3 requires a named othercolor colormap.');
            c=othercolor(clmap,4);
            c(2,1)=0.70;c(2,2)=0.70;c(2,3)=0.70; 
    end
    if neutral_lowest
        % neutral brain colour: grey
        c(1,1)=0.98;c(1,2)=0.98;c(1,3)=0.98;
    end
    colormap(c)
