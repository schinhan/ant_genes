#!/bin/bash

# downsamples and warps the lausanne parcellation (cifti with 32k resolution) to fsaverage5 space with 10k resolution.
# purpose: ABAGEN requires this resolution as mask

wb_command -cifti-separate ../../data/parcellation/Lausanne250.dlabel.nii COLUMN -label CORTEX_LEFT ../../data/parcellation/laus_left.label.gii -label CORTEX_RIGHT ../../data/parcellation/laus_right.label.gii 
wb_command -label-resample ../../data/parcellation/laus_left.label.gii ../../data/templates/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii  ../../data/templates/fsaverage5_std_sphere.L.10k_fsavg_L.surf.gii ADAP_BARY_AREA ../../data/parcellation/Lausanne.250.L.10k_fsavg_5.label.gii -area-metrics ../../data/templates/fs_LR.L.midthickness_va_avg.32k_fs_LR.shape.gii ../../data/templates/fsaverage5.L.midthickness_va_avg.10k_fsavg_L.shape.gii
wb_command -label-resample ../../data/parcellation/laus_right.label.gii ../../data/templates/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii  ../../data/templates/fsaverage5_std_sphere.R.10k_fsavg_R.surf.gii ADAP_BARY_AREA ../../data/parcellation/Lausanne.250.R.10k_fsavg_5.label.gii -area-metrics ../../data/templates/fs_LR.R.midthickness_va_avg.32k_fs_LR.shape.gii ../../data/templates/fsaverage5.R.midthickness_va_avg.10k_fsavg_R.shape.gii