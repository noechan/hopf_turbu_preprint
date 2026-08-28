% combine_ADNI3_export.m
% This script loads ADNI3 .mat files for AD, HC, and MCI groups,
% combines the data, appends a group label, and exports to Excel and CSV.

% Clear workspace
clear;
clc;
sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();
if ~isfolder(P.prepare_metadata), mkdir(P.prepare_metadata); end
cd(P.prepare_metadata)
% Define file names
files = {
    'PTID_ADNI3_AD_MPRAGE_IRFSPGR_all.mat', 'AD';
    'PTID_ADNI3_HC_MPRAGE_IRFSPGR_all.mat', 'HC';
    'PTID_ADNI3_MCI_MPRAGE_IRFSPGR_all.mat', 'MCI'
};

% Collect each group table before concatenating once.
groupTables = cell(size(files, 1), 1);

% Loop over each file
for i = 1:size(files, 1)
    fileName = files{i, 1};
    groupLabel = files{i, 2};
    
    % Load .mat file
    S = load(fileName);
    
    % Extract variable name (assumes first non-meta field is data)
    vars = fieldnames(S);
    data = S.(vars{1});
    
    % Convert to table if it's not already
    if ~istable(data)
        data = array2table(data);
    end
    
    % Append group label
    data.Group = repmat({groupLabel}, height(data), 1);
    
    groupTables{i} = data;
end

combinedData = vertcat(groupTables{:});

% Export to Excel
writetable(combinedData, 'Combined_ADNI3_MPRAGE_IRFSPGR_N238rev.xlsx');

% Export to CSV
writetable(combinedData, 'Combined_ADNI3_MPRAGE_IRFSPGR_N238rev.csv');

disp('Data exported successfully as Excel and CSV files.');
