# ===========================================
# ======= create gwas manhattan plot ========
# ===========================================

# set working directory
cd /slow/projects/coco_genes

# load conda environment
if [ ! -d "envs/default" ]; then
  conda env create --file envs/default.yml -p envs/default
fi
conda activate envs/default

# settings
pFile="results/perm.pvals.txt"
refseqFile="data/refseq/GRCh37_latest_genomic.edit.gff.gz"
targetDir="results/manhattan"
yend=10
ysteps=2
width=11
height=4

# run script
Rscript code/manhattan.R "${pFile}" "${refseqFile}" "${targetDir}" "${yend}" "${ysteps}" "${width}" "${height}"

# ====================================
# ======= create gwas qq-plot ========
# ====================================

# set working directory
cd /slow/projects/coco_genes

# load conda environment
if [ ! -d "envs/default" ]; then
  conda env create --file envs/default.yml -p envs/default
fi
conda activate envs/default

# general settings
pFile='results/perm.pvals.txt' # file with header: gene p.alert p.control p.orient
targetDir="results/qqplot/" # target directory
prune=FALSE # ommit overplotting of variants with pval > 0.01
drawCI=FALSE # draw confidence interval?
xend=6 # upper x axis limit
xsteps=2 # x axis breaks
yend=10 # upper y axis limit
ysteps=2 # y axis breaks
width=3 # plot width (3 inch)
height=4 # plot height (4 inch)

# create plots
Rscript code/qqplot.R "${pFile}" "${targetDir}" "${prune}" "${drawCI}" "${xend}" "${xsteps}" "${yend}" "${ysteps}" "${width}" "${height}"

# =======================================================
# === Combine Manhattan and qq-plots for main article ===
# =======================================================

# set working directory
cd /slow/projects/coco_genes

# load conda environment
if [ ! -d "envs/default" ]; then
  conda env create --file envs/default.yml -p envs/default
fi
conda activate envs/default

# settings
traits="alert,orient,control"
manhattanPlots="results/manhattan/alert.manhattan.png,results/manhattan/orient.manhattan.png,results/manhattan/control.manhattan.png"
qqPlots="results/qqplot/alert.qqplot.png,results/qqplot/orient.qqplot.png,results/qqplot/control.qqplot.png"
outputFile="results/manh.qq.png"
width=14
height=12

# run analysis
Rscript code/manhattan.qq.combine.R "${traits}" "${manhattanPlots}" "${qqPlots}" "${outputFile}" "${width}" "${height}"

