%% coccon_extract_coord_sphere
% gets centroid coordinates of our Lausanne atlas projected on a sphere
% the coordinates are later rotated to create null models

%% requirements:
% templates: fs_L-to-fs_LR_fsaverage sphere, distributed with HCPpipelines
% system(['git clone https://github.com/Washington-University/HCPpipelines.git'])
% templates: HCP's sphere, distributed with HCP_S1200_Atlas_pkXDZ (https://balsa.wustl.edu/reference/pkXDZ, see hcp_datauseterms.txt)
% wb_command
system(['git clone https://github.com/Washington-University/workbench.git'])
% gifti toolbox for matlab
system(['git clone https://github.com/nno/matlab_GIfTI.git'])


%% prepare

% convert Lausanne atlas from cifti to gifti
system(['/Applications/workbench/bin_macosx64/wb_command -cifti-separate ../../data/parcellation/Lausanne250.dlabel.nii COLUMN -label CORTEX_LEFT ../../data/parcellation/Lausanne250.L.label.gii'])
system(['/Applications/workbench/bin_macosx64/wb_command -cifti-separate ../../data/parcellation/Lausanne250.dlabel.nii COLUMN -label CORTEX_RIGHT ../../data/parcellation/Lausanne250.R.label.gii'])

% resample giftis to sphere
system(['/Applications/workbench/bin_macosx64/wb_command -label-resample ../../data/parcellation/Lausanne250.L.label.gii ../../data/templates/S1200.L.sphere.32k_fs_LR.surf.gii ../../data/templates/fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii BARYCENTRIC ../../data/parcellation/Lausanne250.L.fsaverage164.label.gii'])
system(['/Applications/workbench/bin_macosx64/wb_command -label-resample ../../data/parcellation/Lausanne250.R.label.gii ../../data/templates/S1200.R.sphere.32k_fs_LR.surf.gii ../../data/templates/fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.surf.gii BARYCENTRIC ../../data/parcellation/Lausanne250.R.fsaverage164.label.gii'])


%% load giftis (left and right)
lg=gifti('../../data/parcellation/Lausanne250.L.fsaverage164.label.gii');
rg=gifti('../../data/parcellation/Lausanne250.R.fsaverage164.label.gii');

%% load underlying geometry 
lt=gifti('../../data/templates/fs_L-to-fs_LR_fsaverage.L_LR.spherical_std.164k_fs_L.surf.gii'); 
rt=gifti('../../data/templates/fs_R-to-fs_LR_fsaverage.R_LR.spherical_std.164k_fs_R.surf.gii');

%% get mean coordinates

for i=1:length(lg.labels.key)
    slabel(i)=lg.labels.name(i);
    scentroid_l(i,:)=mean(lt.vertices(lg.cdata==lg.labels.key(i),:));
    scentroid_r(i,:)=mean(rt.vertices(rg.cdata==rg.labels.key(i),:));
end

%% clean up
%remove corpus callosum
scentroid_l(contains(slabel,'callo')==1,:)=[];
scentroid_r(contains(slabel,'callo')==1,:)=[];
slabel(contains(slabel,'callo')==1)=[];

%remove ???
scentroid_l(contains(slabel,'?')==1,:)=[];
scentroid_r(contains(slabel,'?')==1,:)=[];
slabel(contains(slabel,'?')==1)=[];

%% re-order to CATO
anttrans=[1	35	2	3	36	37	4	5	6	38	39	40	7	41	42	43	44	8	45	46	47	9	10	48	49	50	51	11	52	53	54	12	55	56	57	13	58	14	59	60	61	15	16	62	17	63	18	19	20	21	64	65	66	67	68	69	22	70	23	71	72	73	74	75	76	77	24	78	79	80	81	25	26	82	83	84	85	86	27	87	88	89	90	91	92	93	94	28	95	96	97	98	99	100	29	101	102	103	104	30	105	106	107	108	31	32	33	34	109	110	111	112	113	114	146	147	115	148	116	117	149	150	151	118	152	153	154	155	156	119	157	158	159	120	121	160	161	162	163	122	164	165	166	123	167	168	124	169	170	125	171	172	173	126	127	174	175	128	176	129	130	177	131	178	132	179	180	181	182	133	183	134	184	185	186	187	188	135	189	190	191	192	136	137	193	194	195	196	197	138	198	199	200	201	202	203	204	139	205	206	207	208	209	210	140	211	212	213	214	141	215	216	217	142	143	144	145	218	219];

slabel=slabel(anttrans);
scentroid_l=scentroid_l(anttrans,:);
scentroid_r=scentroid_r(anttrans,:);

%% remove other hemisphere
scentroid_l(isnan(mean(scentroid_l,2)),:)=[];
scentroid_r(isnan(mean(scentroid_r,2)),:)=[];

%% save results
label_left=slabel(contains(slabel,'L'));
label_right=slabel(contains(slabel,'R'));
coord_left=scentroid_l;
coord_right=scentroid_r;

savefile='../../data/parcellation/sphere_coordinates.mat';
save(savefile,'coord_right','coord_left',"label_right","label_left")

