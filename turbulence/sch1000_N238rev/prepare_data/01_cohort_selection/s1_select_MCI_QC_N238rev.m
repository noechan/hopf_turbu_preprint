clearvars
sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();
% MATLAB script to filter subjects and save indices by batch
MCI_data_path = P.participants_cross_sectional;
code_path = P.prepare_metadata;
if ~isfolder(code_path), mkdir(code_path); end
% Specify the file name
filename = 'ADNI3_MCI_ALL_MRI_sc_analysis_N238rev_idx4MRIselection.xlsx';
sheet = 'MCI_removed_N238rev_idx'; % Specify the sheet name
% Read the data from the Excel file
data = readtable(fullfile(MCI_data_path,filename), 'Sheet', sheet);

% Initialize containers for connIDs variables by batch
unique_batches = unique(data.CONNBatch); % Extract unique batch names
connIDs = struct(); % Use a structure to store variables dynamically
ptid = struct(); % Structure to store PTIDs

% Loop through each batch to filter data and extract connIDs
for i = 1:length(unique_batches)
    batch_name = unique_batches{i}; % Current batch name
    
    % Filter data for the current batch and exclude rows where Total_removed == 1
    batch_data = data(strcmp(data.CONNBatch, batch_name) & data.Total_removed ~= 1, :);
    
    % % Extract the connIDs and PTID
    connIDs_field = ['connIDs_', strrep(batch_name, ' ', '_'),'_MCI']; % Replace spaces with underscores
    ptid_field = ['PTID_', strrep(batch_name, ' ', '_'),'_MCI']; % Replace spaces with underscores
    connIDs.(connIDs_field) = batch_data.CONNID;
    ptid.(ptid_field) = batch_data.PTID; % Extract the PTID column for the batch
end

% Save each batch's connIDs and PTID as separate variables in the workspace
field_names_connIDs = fieldnames(connIDs);
field_names_ptid = fieldnames(ptid);

cd(code_path)

% Save connIDs variables
for i = 1:length(field_names_connIDs)
    varname = field_names_connIDs{i};
    assignin('base', varname, connIDs.(varname)); % Assign to the base workspace
    save([varname, '.mat'], varname); % Save as a .mat file
end

% Save PTID variables
for i = 1:length(field_names_ptid)
    varname = field_names_ptid{i};
    assignin('base', varname, ptid.(varname)); % Assign to the base workspace
    save([varname, '.mat'], varname); % Save as a .mat file
end
