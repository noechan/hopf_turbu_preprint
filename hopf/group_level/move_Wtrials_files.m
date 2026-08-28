clear; clc;

% Define the parent directory containing all WG files and group folders
parentDir = '/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/Wtrials/';  
% Define group subfolders
groupFolders = {'HC_ABneg', 'HC_ABpos', 'MCI_ABpos', 'AD_ABpos'};

% Create subfolders if they do not exist
for i = 1:length(groupFolders)
    folderPath = fullfile(parentDir, groupFolders{i});
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
        fprintf('Created folder: %s\n', folderPath);
    end
end

% Get all WG files
wgFiles = dir(fullfile(parentDir, 'Wtrials_*.mat'));

% Loop through each file and move based on the group label
for i = 1:length(wgFiles)
    fileName = wgFiles(i).name;
    srcPath  = fullfile(parentDir, fileName);

    % Determine destination folder based on pattern in filename
    if contains(fileName, 'HC_ABneg')
        destFolder = 'HC_ABneg';
    elseif contains(fileName, 'HC_ABpos')
        destFolder = 'HC_ABpos';
    elseif contains(fileName, 'MCI_ABpos')
        destFolder = 'MCI_ABpos';
    elseif contains(fileName, 'AD_ABpos')
        destFolder = 'AD_ABpos';
    else
        fprintf('Skipping file (no match found): %s\n', fileName);
        continue; % skip if no match
    end

    % Move file
    destPath = fullfile(parentDir, destFolder, fileName);
    movefile(srcPath, destPath);
    fprintf('Moved %s → %s\n', fileName, destFolder);
end

fprintf('\n✅ File transfer completed successfully.\n');
