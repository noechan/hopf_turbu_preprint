
clear; clc
code_path='/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/';
outdir_HC_ABneg='/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/G_range/HC_ABneg/';
outdir_HC_ABpos='/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/G_range/HC_ABpos/';
outdir_MCI_ABpos='/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/G_range/MCI_ABpos/';
outdir_AD_ABpos='/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/G_range/AD_ABpos/';
addpath(outdir_HC_ABneg); addpath(outdir_HC_ABpos);addpath(outdir_MCI_ABpos);addpath(outdir_AD_ABpos)
addpath(genpath(code_path))
G=100;
G_range=0.:0.01:3;

cond=1;
cond_HC_ABneg=load_condition_HC_ABneg(G,cond);
[minerrhete_HC_ABneg,ioptG_HC_ABneg]=min(mean(cond_HC_ABneg.err_hete_range,2));
optG_HC_ABneg=G_range(ioptG_HC_ABneg);

cd(fullfile(outdir_HC_ABneg,'optG_HC_ABneg'))
save('optG_HC_ABneg.mat','optG_HC_ABneg')

cond=1;
cond_HC_ABpos=load_condition_HC_ABpos(G,cond);
[minerrhete_HC_ABpos,ioptG_HC_ABpos]=min(mean(cond_HC_ABpos.err_hete_range,2));
optG_HC_ABpos=G_range(ioptG_HC_ABpos);

cd(fullfile(outdir_HC_ABpos,'optG_HC_ABpos'))
save('optG_HC_ABpos.mat','optG_HC_ABpos')

cond=1;
cond_MCI_ABpos=load_condition_MCI_ABpos(G,cond);
[minerrhete_MCI_ABpos,ioptG_MCI_ABpos]=min(mean(cond_MCI_ABpos.err_hete_range,2));
optG_MCI_ABpos=G_range(ioptG_MCI_ABpos);

cd(fullfile(outdir_MCI_ABpos,'optG_MCI_ABpos'))
save('optG_MCI_ABpos.mat','optG_MCI_ABpos')

cond=1;
cond_AD_ABpos=load_condition_AD_ABpos(G,cond);
[minerrhete_AD_ABpos,ioptG_AD_ABpos]=min(mean(cond_AD_ABpos.err_hete_range,2));
optG_AD_ABpos=G_range(ioptG_AD_ABpos);

cd(fullfile(outdir_AD_ABpos,'optG_AD_ABpos'))
save('optG_AD_ABpos.mat','optG_AD_ABpos')

% Plot selected G range for comparative purposes
x=G_range(1:50)';
figure;
shadedErrorBar(x, cond_HC_ABneg.err_hete_range(1:50,:)',{@mean,@std},'lineprops',{'k-','color',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7]}); %grey
hold on
shadedErrorBar(x, cond_HC_ABpos.err_hete_range(1:50,:)',{@mean,@std},'lineprops',{'k-','color',[0.4940 0.1840 0.5560],'markerfacecolor',[0.4940 0.1840 0.5560]}); %purple
shadedErrorBar(x, cond_MCI_ABpos.err_hete_range(1:50,:)',{@mean,@std},'lineprops',{'k-','color',[0.8500 0.3250 0.0980],'markerfacecolor',[0.8500 0.3250 0.0980]}); %orange
shadedErrorBar(x, cond_AD_ABpos.err_hete_range(1:50,:)',{@mean,@std},'lineprops',{'k-','color',[0 0.4470 0.7410],'markerfacecolor',[0 0.4470 0.7410]}); %blue
ylabel('Fitting');xlabel('Global coupling G')
xline(optG_HC_ABneg, 'Color',[0.7 0.7 0.7], 'LineWidth',1.5)
xline(optG_HC_ABpos, 'Color',[0.4940 0.1840 0.5560], 'LineWidth',1.5)
xline(optG_MCI_ABpos, 'Color',[0.8500 0.3250 0.0980], 'LineWidth',1.5)
xline(optG_AD_ABpos, 'Color',[0 0.4470 0.7410], 'LineWidth',1.5)
legend({'HC_ABneg','HC_ABpos', 'MCI_ABpos', 'AD_ABpos'}, 'Location','northeastoutside')
xticks([0 0.1 0.2 0.3 0.4 0.5])


cd('/path/to/ADNI3/code/HPC_Hopf_DTI_1000_Staging/G_range/')
print('Global coupling Fitting ADNI3 Staging 4(N=238)','-dpng')
saveas(gcf,'Global coupling Fitting ADNI3 Staging 4 (N=238)','fig')