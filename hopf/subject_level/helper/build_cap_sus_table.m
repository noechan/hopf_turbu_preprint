function T = build_cap_sus_table(pertFile, empFile)
% BUILD_CAP_SUS_TABLE
% Creates a table with PTID, CAP and Susceptibility
% using two different .mat files.

%% ---------------- Load files ----------------
Spert = load(pertFile);
Semp  = load(empFile);

%% ---------------- Find PTID (from empirical file) ----------------
fn = fieldnames(Semp);

PTID = [];
for i = 1:numel(fn)
    if contains(fn{i},'PTID','IgnoreCase',true)
        PTID = Semp.(fn{i});
        fprintf('PTID found in empirical file: %s\n',fn{i});
        break
    end
end

if isempty(PTID)
    error('No PTID found in empirical correlation file.')
end

PTID = normalize_ptid(PTID);
N = numel(PTID);

%% ---------------- Find CAP & Susceptibility (perturbation file) ----------------

[CAP,capName] = find_vector(Spert,["infocap"],N);
[SUS,susName] = find_vector(Spert,["suscep"],N);

fprintf('CAP variable: %s\n',capName);
fprintf('Susceptibility variable: %s\n',susName);

%% ---------------- Alignment safety check ----------------
% (this is CRITICAL for ADNI)

% Try to find PTID also inside perturbation file
pertPTID = [];
fn2 = fieldnames(Spert);
for i = 1:numel(fn2)
    if contains(fn2{i},'PTID','IgnoreCase',true)
        pertPTID = normalize_ptid(Spert.(fn2{i}));
        break
    end
end

if ~isempty(pertPTID)

    % reorder CAP/SUS to match empirical PTID order
    [commonID, idxEmp, idxPert] = intersect(PTID,pertPTID,'stable');

    if length(commonID) ~= N
        error('Subjects mismatch between files.');
    end

    CAP = CAP(idxPert);
    SUS = SUS(idxPert);
end

%% ---------------- Build table ----------------

T = table(PTID, CAP(:), SUS(:), ...
          'VariableNames',{'PTID','CAP','Susceptibility'});

end


%% ================== helper functions ==================

function [saidVec,varName] = find_vector(S,patterns,N)

fn = fieldnames(S);
saidVec = [];
varName = "";

for i = 1:numel(fn)

    name = string(fn{i});
    val  = S.(fn{i});

    for p = 1:numel(patterns)
        if contains(lower(name),lower(patterns(p)))

            if isnumeric(val) && numel(val)==N
                saidVec = val(:);
                varName = name;
                return
            end

        end
    end
end

error('Could not find variable matching patterns: %s',strjoin(patterns,", "));
end


function PTID = normalize_ptid(raw)

if isstring(raw)
    PTID = raw(:);
elseif iscell(raw)
    PTID = string(raw(:));
elseif ischar(raw)
    PTID = string(cellstr(raw));
else
    error('Unsupported PTID format')
end

end

