%%  ANT second level ciftis
% This code takes the group level activation maps from the attention
% network test (alerting, validity/orienting, control), converts the maps
% to gifti, loads the lausanne parcellation as gifti, and overlays the two
% data to derive mean activation estimates for each parcel. 
%% requirements
% wb_command
system(['git clone https://github.com/Washington-University/workbench.git'])
% gifti toolbox for matlab
system(['git clone https://github.com/nno/matlab_GIfTI.git'])
%%
% get file names
x=dir('../../data/fMRI/*zTstat*');

% separate cifti to gifti (using wb_command)
system(['/Applications/workbench/bin_macosx64/wb_command -cifti-separate ' [x(1).folder '/' x(1).name] ' COLUMN  -metric CORTEX_LEFT ' x(1).folder '/alert.left.func.gii']); 
system(['/Applications/workbench/bin_macosx64/wb_command -cifti-separate ' [x(1).folder '/' x(1).name] ' COLUMN  -metric CORTEX_RIGHT ' x(1).folder '/alert.right.func.gii']);
system(['/Applications/workbench/bin_macosx64/wb_command -cifti-separate ' [x(4).folder '/' x(4).name] ' COLUMN  -metric CORTEX_LEFT ' x(4).folder '/orient.left.func.gii']); 
system(['/Applications/workbench/bin_macosx64/wb_command -cifti-separate ' [x(4).folder '/' x(4).name] ' COLUMN  -metric CORTEX_RIGHT ' x(4).folder '/orient.right.func.gii']);
system(['/Applications/workbench/bin_macosx64/wb_command -cifti-separate ' [x(5).folder '/' x(5).name] ' COLUMN  -metric CORTEX_LEFT ' x(5).folder '/control.left.func.gii']); 
system(['/Applications/workbench/bin_macosx64/wb_command -cifti-separate ' [x(5).folder '/' x(5).name] ' COLUMN  -metric CORTEX_RIGHT ' x(5).folder '/control.right.func.gii']);

% list file names
y=dir('../../data/fMRI/*gii');

% loaad giftis
for i=1:3
    g=gifti([y(i).folder '/' y(i).name]);
    antleft(:,i)=g.cdata;
    g=gifti([y(i+2).folder '/' y(i).name]);
    antright(:,i)=g.cdata;
end

%% load Laus-250 parcellation (as gifti)
L=gifti('../../data/parcellation/left_laus250.label.gii')
R=gifti('../../data/parcellation/right_laus250.label.gii');
laus(:,1)=L.cdata; laus(:,2)=R.cdata; 

%% extract mean activation from parcels
% unique labels, without trailing zero
ul=unique(laus(:,1)); ur=unique(laus(:,2)); ul(1)=[];ur(1)=[];
for j=1:3
    for i=1:length(ul)
        leftvec(i,j)=mean(antleft(laus(:,1)==ul(i),j));
    end
    for i=1:length(ur)
        rightvec(i,j)=mean(antright(laus(:,2)==ur(i),j));
    end
end

%% pull anatomical labels
for i=1:length(ul)
    lleft(i)=L.labels.name(L.labels.key==(ul(i)));
end
for i=1:length(ur)
    lright(i)=R.labels.name(R.labels.key==(ur(i)));
end

%% combine both hemispheres
ant_lausanne=vertcat(leftvec,rightvec);
labels=horzcat(lleft,lright);

%% save as .mat and .csv
savefile='../../data/fMRI/ant_lausanne.mat';
save(savefile,'ant_lausanne','labels')

alert=ant_lausanne(:,1); orient=ant_lausanne(:,2); control=ant_lausanne(:,3);
T=table(labels', alert,orient,control);
writetable(T,'../../data/fMRI/ant_lausanne.csv');
