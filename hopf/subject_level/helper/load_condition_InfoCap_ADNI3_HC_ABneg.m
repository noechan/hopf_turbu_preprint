function [vcondition]=load_condition_InfoCap_ADNI3_HC_ABneg(trial,cond)

for i=1:trial
    load(sprintf('Wtrials_%03d_%d_err_hete_ADNI3_HC_ABneg_GEC_sch1000.mat',i,cond));
    vcondition.infocap_all(i,:)=infocapacity;
    vcondition.suscep_all(i,:)=susceptibility;
end