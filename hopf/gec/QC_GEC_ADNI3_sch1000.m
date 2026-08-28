%% Quality Control

% 1) Global EC asymmetry
C = Ceffgroup;
Ag = norm((C - C.')/2,'fro') / max(norm(C,'fro'), eps);
fprintf('Global EC asymmetry: %.4f\n', Ag);


% 2) Residual asymmetry at τ=1 (should be < ~0.178 baseline; ideally ~0.12–0.15)
Rplus = (COVtauemp - COVtausim);
Ranti = (Rplus - Rplus.')/2;
Ag_res = norm(Ranti,'fro') / max(norm(Rplus,'fro'), eps);
fprintf('Residual asymmetry (τ=1): %.4f\n', Ag_res);

%Residual asymmetry (τ=1): 0.1768

% 3) FC & lag(1TR) fits (expect FC r ~0.6–0.8; R² ~0.3–0.5; Lag r ~0.5–0.7; R² ~0.2–0.4)
N = size(Ceffgroup,1); m = triu(true(N),1);
fc_emp = FCemp(m); fc_sim = FCsim(m);
lag_emp = COVtauemp(m); lag_sim = COVtausim(m);
r_fc  = corr(fc_emp, fc_sim);   R2_fc  = 1 - sum((fc_emp-fc_sim).^2)/sum((fc_emp-mean(fc_emp)).^2);
r_lag = corr(lag_emp, lag_sim); R2_lag = 1 - sum((lag_emp-lag_sim).^2)/sum((lag_emp-mean(lag_emp)).^2);
fprintf('FC: r=%.3f R^2=%.3f | Lag(1TR): r=%.3f R^2=%.3f\n', r_fc, R2_fc, r_lag, R2_lag);

%HC: FC: r=0.689 R^2=0.349 | Lag(1TR): r=0.654 R^2=0.317

% 4) SC–EC correlation (healthy band ~0.4–0.6)
SC = sc_schaefer; r_sc = corr(SC(m), Ceffgroup(m));
fprintf('SC–EC correlation: r=%.3f\n', r_sc);

%SC–EC correlation: r=0.776

% 5) Stability (must be negative)
lam = eig(A); fprintf('max Re(lambda(A))=%.3f\n', max(real(lam)));

% 6)  SC–FC correlation ---
N = size(FCemp,1);                  % number of parcels
maskUT = triu(true(N),1);           % upper triangle without diagonal

% Vectorize upper triangle of SC and FC
sc_vec = sc_schaefer_400(maskUT);
fc_vec = FCemp(maskUT);

% Correlation
r_sc_fc = corr(sc_vec, fc_vec);

fprintf('SC–FC correlation: r = %.3f\n', r_sc_fc);

% 7) IMPORTANT!! Anti-symmetric alignment check
N = size(Ceffgroup,1);

% Anti-symmetric part of EC
AEC = (Ceffgroup - Ceffgroup.')/2;

% Anti-symmetric part of empirical lag covariance (τ = 1)
ALag = (COVtauemp - COVtauemp.')/2;

% Upper triangle mask (exclude diagonal)
maskUT = triu(true(N),1);

% Vectorize and correlate
r_align = corr(AEC(maskUT), ALag(maskUT));

fprintf('Alignment EC asymmetry vs lag asymmetry (τ=1): r = %.3f\n', r_align);
