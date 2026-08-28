%% Turbulence observables: cognitively healthy Aβ-positive group
clearvars; close all

sch1000_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(sch1000_root, '-begin');
P = setup_sch1000_paths();

input_file = fullfile(P.input.HC, ...
    'tseries_ADNI3_HC_ABetaPos_MPRAGE_IRFSPGR_sch1000_N238rev.mat');
input_data = load(input_file, 'tseries_HC_ABpos', 'PTID_HC_ABpos');
xs = input_data.tseries_HC_ABpos(:);
PTID = cellstr(string(input_data.PTID_HC_ABpos(:)));
group_key = 'HC_ABetaPos';
group_label = 'HC Aβ+';
output_dir = P.output.HC_ABetaPos;

NPARCELS = 1000;
Tmax = 197;
if numel(xs) ~= numel(PTID) || isempty(xs)
    error('%s input has %d time series and %d PTIDs.', ...
        group_label, numel(xs), numel(PTID));
end
for sub = 1:numel(xs)
    if ~isnumeric(xs{sub}) || ~isequal(size(xs{sub}), [NPARCELS Tmax])
        error('%s participant %d (%s) has size %s; expected %d x %d.', ...
            group_label, sub, PTID{sub}, mat2str(size(xs{sub})), ...
            NPARCELS, Tmax);
    end
end

load(fullfile(P.sch1000_root, 'SchaeferCOG.mat'), 'SchaeferCOG');
NSUB = numel(xs);
TR = 3;
NR = 400;
NRini = 20;
NRfin = 80;
LAMBDA = [0.27 0.24 0.21 0.18 0.15 0.12 0.09 0.06 0.03 0.01];
NLAMBDA = numel(LAMBDA);

fprintf('Calculating %s: %d participants.\n', group_label, NSUB);
fprintf('Input: %s\n', input_file);

fnq = 1 / (2 * TR);
Wn = [0.008 / fnq, 0.08 / fnq];
[bfilt, afilt] = butter(2, Wn);

rr = zeros(NPARCELS, NPARCELS);
for i = 1:NPARCELS
    for j = 1:NPARCELS
        rr(i,j) = norm(SchaeferCOG(i,:) - SchaeferCOG(j,:));
    end
end
delta = max(rr, [], 'all') / NR;
xrange = delta / 2 + delta * (0:NR-1);

C1 = zeros(NLAMBDA, NPARCELS, NPARCELS);
for ilam = 1:NLAMBDA
    C1(ilam,:,:) = exp(-LAMBDA(ilam) * rr);
end

Turbulence_global_sub = zeros(NLAMBDA, NSUB);
Turbulence_node_sub = zeros(NLAMBDA, NPARCELS, NSUB);
Local_KoP_sub = zeros(NLAMBDA, NSUB);
TransferLambda_sub = zeros(NLAMBDA, NSUB);
InformationCascade_sub = zeros(1, NSUB);
Transfer_sub = zeros(NLAMBDA, NSUB);
gKoP = zeros(1, NSUB);
Meta = zeros(1, NSUB);

signal_filt = zeros(NPARCELS, Tmax);
Phases = zeros(NPARCELS, Tmax);
LocalOrderParameter = zeros(NLAMBDA, NPARCELS, Tmax);

for sub = 1:NSUB
    fprintf('Participant %d/%d: %s\n', sub, NSUB, PTID{sub});
    ts = xs{sub};

    for seed = 1:NPARCELS
        seed_ts = detrend(ts(seed,:) - mean(ts(seed,:)));
        signal_filt(seed,:) = filtfilt(bfilt, afilt, seed_ts);
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end

    complex_phase = complex(cos(Phases), sin(Phases));
    for ilam = 1:NLAMBDA
        C1lam = squeeze(C1(ilam,:,:));
        for i = 1:NPARCELS
            weighted_phase = C1lam(i,:)' .* complex_phase;
            local_order = abs(sum(weighted_phase, 1, 'omitnan') / ...
                sum(C1lam(i,:)));
            LocalOrderParameter(ilam,i,:) = local_order;
            Turbulence_node_sub(ilam,i,sub) = std(local_order, 0, 'omitnan');
        end
        all_local_order = reshape(LocalOrderParameter(ilam,:,:), [], 1);
        Turbulence_global_sub(ilam,sub) = std(all_local_order, 0, 'omitnan');
        Local_KoP_sub(ilam,sub) = mean(all_local_order, 'omitnan');
    end

    global_order = abs(sum(complex_phase, 1)) / NPARCELS;
    gKoP(sub) = mean(global_order, 'omitnan');
    Meta(sub) = std(global_order, 0, 'omitnan');

    for ilam = 1:NLAMBDA-1
        current_scale = squeeze(LocalOrderParameter(ilam,:,1:end-1))';
        next_scale = squeeze(LocalOrderParameter(ilam+1,:,2:end))';
        [cc, pp] = corr(next_scale, current_scale);
        TransferLambda_sub(ilam+1,sub) = mean( ...
            abs(cc(pp < 0.05)), 'omitnan');
    end
    InformationCascade_sub(sub) = mean( ...
        TransferLambda_sub(2:NLAMBDA,sub), 'omitnan');

    for ilam = 1:NLAMBDA
        fclam = corrcoef(squeeze(LocalOrderParameter(ilam,:,:))');
        numind = zeros(1, NR);
        fcra = zeros(1, NR);
        for i = 1:NPARCELS
            for j = 1:NPARCELS
                if i == j
                    continue
                end
                index = floor(rr(i,j) / delta) + 1;
                if index == NR + 1
                    index = NR;
                end
                if ~isnan(fclam(i,j))
                    fcra(index) = fcra(index) + fclam(i,j);
                    numind(index) = numind(index) + 1;
                end
            end
        end

        grandcorrfcn = fcra ./ numind;
        available_bins = find(~isnan(grandcorrfcn));
        fit_bins = available_bins(available_bins > NRini & available_bins < NRfin);
        xcoor = zeros(1, numel(fit_bins));
        ycoor = zeros(1, numel(fit_bins));
        nn = 0;
        for kk = 1:numel(fit_bins)
            radial_bin = fit_bins(kk);
            if grandcorrfcn(radial_bin) > 0
                nn = nn + 1;
                xcoor(nn) = log(xrange(radial_bin));
                ycoor(nn) = log(grandcorrfcn(radial_bin) / ...
                    grandcorrfcn(fit_bins(1)));
            end
        end
        xcoor = xcoor(1:nn);
        ycoor = ycoor(1:nn);

        linfunc = @(A, x) A(1) * x + A(2);
        options = optimset('MaxFunEvals', 10000, 'MaxIter', 1000, ...
            'Display', 'final');
        Afit = lsqcurvefit(linfunc, [-1 1], xcoor, ycoor, ...
            [-4 -10], [4 10], options);
        Transfer_sub(ilam,sub) = abs(Afit(1));
    end
end

if ~isfolder(output_dir)
    mkdir(output_dir);
end
save(fullfile(output_dir, ...
    'turbu_all_measurements_ADNI3_HC_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000.mat'), ...
    'group_key', 'input_file', 'PTID', 'LAMBDA', ...
    'Turbulence_global_sub', 'Local_KoP_sub', 'Transfer_sub', ...
    'InformationCascade_sub', 'TransferLambda_sub', 'gKoP', 'Meta', '-v7.3');
save(fullfile(output_dir, ...
    'turbu_by_node_ADNI3_HC_ABetaPos_MPRAGE_IRFSPGR_N238rev_sch1000.mat'), ...
    'group_key', 'input_file', 'PTID', 'Turbulence_node_sub', '-v7.3');

fprintf('Finished %s. Results saved in %s\n', group_label, output_dir);
