clear all
addpath(genpath('/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/'));
addpath('/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/'); %codepath
outdir='/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/plots/';
trials=100;
groups={'HC_ABneg','HC_ABpos','MCI_ABpos','AD_ABpos'};

HC_ABneg=load_condition_InfoCap_ADNI3_HC_ABneg(trials,1);
HC_ABpos=load_condition_InfoCap_ADNI3_HC_ABpos(trials,1);
MCI_ABpos=load_condition_InfoCap_ADNI3_MCI_ABpos(trials,1);
AD_ABpos=load_condition_InfoCap_ADNI3_AD_ABpos(trials,1);

%% To export for further stats in R 
HC_ABneg_infocap=HC_ABneg.infocap_all;
HC_ABpos_infocap=HC_ABpos.infocap_all;
MCI_ABpos_infocap=MCI_ABpos.infocap_all;
AD_ABpos_infocap=AD_ABpos.infocap_all;

HC_ABneg_suscep=HC_ABneg.suscep_all;
HC_ABpos_suscep=HC_ABpos.suscep_all;
MCI_ABpos_suscep=MCI_ABpos.suscep_all;
AD_ABpos_suscep=AD_ABpos.suscep_all;

cd(outdir)
save ('perturbation_HC_ABneg.mat', 'HC_ABneg_infocap','HC_ABneg_suscep')
save ('perturbation_HC_ABpos.mat', 'HC_ABpos_infocap','HC_ABpos_suscep')
save ('perturbation_MCI_ABpos.mat', 'MCI_ABpos_infocap','MCI_ABpos_suscep')
save ('perturbation_AD_ABpos.mat', 'AD_ABpos_infocap','AD_ABpos_suscep')


%% Visualization
figure
boxplot([HC_ABneg.infocap_all HC_ABpos.infocap_all MCI_ABpos.infocap_all AD_ABpos.infocap_all],'labels',groups)
title('Information capacity');% remove NaNs just in case

[pval, observeddifference, effectsize]  = permutationTest((HC_ABneg.infocap_all),(HC_ABpos.infocap_all), 1000); 
p1P(1,1)=pval; clear pval
[pval, observeddifference, effectsize]  = permutationTest((HC_ABpos.infocap_all),(MCI_ABpos.infocap_all), 1000); 
p1P(1,2)=pval; clear pval
[pval, observeddifference, effectsize]  = permutationTest((MCI_ABpos.infocap_all),(AD_ABpos.infocap_all), 1000); 
p1P(1,3)=pval; clear pval
[pval, observeddifference, effectsize]  = permutationTest((HC_ABneg.infocap_all),(MCI_ABpos.infocap_all), 1000); 
p1P(1,4)=pval; clear pval
[pval, observeddifference, effectsize]  = permutationTest((HC_ABneg.infocap_all),(AD_ABpos.infocap_all), 1000); 
p1P(1,5)=pval; clear pval
H=sigstar({[1,2],[2,3],[3,4],[1,3],[1,4]},p1P);

print('Info Capacity ADNI3 N=238 Permutation 1000 both','-dpng')
saveas(gcf,'Info Capacity ADNI3 N=238 Permutation 1000 both','fig')
close()
% FDR correction
[h_t,~,~,adj_p_t] = fdr_bh(p1P(1,1:5),0.05,'pdep','yes'); 
save('InfoCap_rejectedH0.mat','h_t','adj_p_t')

figure
boxplot([HC_ABneg.suscep_all HC_ABpos.suscep_all MCI_ABpos.suscep_all AD_ABpos.suscep_all],'labels',groups)
title('Susceptibility');

[pval, observeddifference, effectsize]  = permutationTest((HC_ABneg.suscep_all),(HC_ABpos.suscep_all), 1000); 
p2P(1,1)=pval; clear pval
[pval, observeddifference, effectsize]  = permutationTest((HC_ABpos.suscep_all),(MCI_ABpos.suscep_all), 1000); 
p2P(1,2)=pval; clear pval
[pval, observeddifference, effectsize]  = permutationTest((MCI_ABpos.suscep_all),(AD_ABpos.suscep_all), 1000); 
p2P(1,3)=pval; clear pval
[pval, observeddifference, effectsize]  = permutationTest((HC_ABneg.suscep_all),(MCI_ABpos.suscep_all), 1000); 
p2P(1,4)=pval; clear pval
[pval, observeddifference, effectsize]  = permutationTest((HC_ABneg.suscep_all),(AD_ABpos.suscep_all), 1000); 
p2P(1,5)=pval; clear pval
H=sigstar({[1,2],[2,3],[3,4],[1,3],[1,4]},p2P);

print('Susceptibility ADNI3 N=238 Permutation 1000 both','-dpng')
saveas(gcf,'Susceptibility ADNI3 N=238 Permutation 1000 both','fig')
close()

% FDR correction
[h_t,~,~,adj_p_t] = fdr_bh(p2P(1,1:5),0.05,'pdep','yes'); % HC- vs HC+
save('InfoSus_rejectedH0.mat','h_t','adj_p_t')

