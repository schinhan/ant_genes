%% code to correlate ANT activation (fMRI) to receptor/transporter availabbility (PET)
% required: Data from Justine Hansen
 system(['git clone https://github.com/netneurolab/hansen_receptors.git data/pet'])
 % also required: rotate_parcellation (by Frantisek Vasa)
 system(['git clone  https://github.com/frantisekvasa/rotate_parcellation.git'])
%% load PET maps, Lausanne 219
% read receptormaps
recep=load('data/pet/hansen_receptors-main/results/receptor_data_scale125.csv');

% load ANT, Lausanne 219
load data/fmri/ant_lausanne.mat

% load transforms
load code/PET/anttrans.mat
load code/PET/pet2cato.mat

%reorder: we reorder all brain maps to match Cato's lausanne template
recep1=recep(firstrans,:);ant_recep=recep1(secondtrans,:); % first trans removes subcortical, second trans aligns with CATO
ant_cons=ant_lausanne(anttrans,:);

% receptor / transporter labels
rlabels={'5HT1a' '5HT1b' '5HT2a' '5HT4' '5HT6' '5HTT' 'A4B2' 'CB1' 'D1' 'D2' 'DAT' 'GABAa' 'H3' 'M1' 'mGluR5' 'MOR' 'NET' 'NMDA' 'VAChT'};

%% compute correlations

for i=1:3
    for j=1:19
        [r_ant_rec(j,i),p_ant_rec(j,i)]=corr(ant_cons(:,i),ant_recep(:,j),'type','pearson');
    end
end


%% permutations

% load centroids from Lausanne250 group parcellation
load data/parcellation/sphere_coordinates.mat

% rotate
for i=1:3
    for j=1:19
        perm_id(:,:,i,j) = rotate_parcellation(coord_left,coord_right,5000);
    end
end

% run permutation tests
targets=[7 9 10 11 17 19]; % those are the ones we have an hypothesis on

for i=1:3
    for j=1:19
        for k=1:5000
            null_corr_all(i,j,k)=corr(ant_cons(perm_id(:,k,i,j),i),ant_recep(:,j),'type','pearson');
        end
    end
end

for i=1:3
    for j=1:19
        pnull_all(j,i,1)=1-(sum(r_ant_rec(j,i)>null_corr_all(i,j,:))/5000);
        pnull_all(j,i,2)=1-(sum(r_ant_rec(j,i)<null_corr_all(i,j,:))/5000);
    end
end

% fdr-correct:
for i=1:3; for j=1:2
    [h(:,i,j),~,pfdr(:,i,j)]=fdr_bh(pnull_all(setdiff(1:19,targets),i,j),.05,'pdep','false');
end; end

%% save results
petmap=rlabels(targets)';
alert_r=r_ant_rec(targets,1); alert_p=pnull_all(targets,1,1);    sys=pnull_all(targets,1,2); alert_p(find(alert_r<0))=sys(find(alert_r<0));
orient_r=r_ant_rec(targets,2); orient_p=pnull_all(targets,2,1);  sys=pnull_all(targets,2,2); orient_p(find(orient_r<0))=sys(find(orient_r<0));
control_r=r_ant_rec(targets,3); control_p=pnull_all(targets,3,1);sys=pnull_all(targets,3,2); control_p(find(control_r<0))=sys(find(control_r<0));
ant_results_hypothesis=table(petmap,alert_r,alert_p,orient_r,orient_p,control_r,control_p)

petmap=rlabels(setdiff(1:19,targets))';
alert_r=r_ant_rec(setdiff(1:19,targets),1); alert_p=pfdr(:,1,1); alert_p(find(alert_r<0))=pfdr(find(alert_r<0),1,2);
orient_r=r_ant_rec(setdiff(1:19,targets),2); orient_p=pfdr(:,2,1); orient_p(find(orient_r<0))=pfdr(find(orient_r<0),2,2);
control_r=r_ant_rec(setdiff(1:19,targets),3); control_p=pfdr(:,3,1); control_p(find(control_r<0))=pfdr(find(control_r<0),3,2);
ant_results_exploratory=table(petmap,alert_r,alert_p,orient_r,orient_p,control_r,control_p)

% for table
savefile='data/pet/ant_pet_results.mat'
save(savefile,'ant_results_exploratory','ant_results_hypothesis')

% for plotting
savefile='visualization/pet/stuff4plot.mat'
save(savefile,'ant_cons','ant_recep','rlabels','targets')

