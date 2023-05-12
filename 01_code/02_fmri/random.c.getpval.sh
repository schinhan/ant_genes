#!/bin/bash

# =============================================
# === get p-values from permutation testing ===
# =============================================

# set working directory
cd /slow/projects/coco_genes

# settings
nworkers=50
npermutations=1000000

# count |expected correlations| larger than |observed correlation|
for cond in alert orient control; do
	echo "Starting with ${cond}"
(
	for ((i=1; i<=$nworkers;i++)); do (
	awk -F'\t' '
		function abs(v) {return v < 0 ? -v : v}
		FNR==1 { next }
		NR==FNR { corr[NR-1]=abs($2); sum[NR-1]=0; geneCount++; next }
		{for (i=1;i<=NF;i++) { if(abs($i)>=corr[i]) { sum[i]++}}}
		END {for (i=1;i<=geneCount;i++) printf "%d%s", sum[i], (i==geneCount?"\n":" ")}
		' results/${cond}.obs.corr.txt results/${cond}/exp.corr.$(printf '%04d' $i).txt > results/${cond}/sum.twotailed.$(printf '%04d' $i).txt
	) &
	done
	wait
)
cat results/${cond}/sum.twotailed*.txt | awk -vnpermutations="${npermutations}" 'NR==1 {geneCount=NF} {for (i=1;i<=NF;i++) {sum[i]+=$i}} END {for (i=1;i<=geneCount;i++) printf "%.6f%s", sum[i]/npermutations, (i==geneCount?"\n":" ")}' > results/${cond}.pval.twotailed.txt
chmod 770 results/${cond}.pval.twotailed.txt
done
