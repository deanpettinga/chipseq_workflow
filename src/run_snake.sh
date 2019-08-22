#PBS -l walltime=100:00:00
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=1
#PBS -M dean.pettinga@vai.org
#PBS -m ae
#PBS -N snake_ChIP
#PBS -o logs/snake_ChIP.o
#PBS -e logs/snake_ChIP.e

cd ${PBS_O_WORKDIR}
# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
snakemake --use-conda -n > runs/dag_${TIME}.txt
snakemake --dag | dot -Tpng > runs/dag_${TIME}.png

snakemake \
-s Snakefile \
-j 48 \
--cluster-config src/cluster.yaml \
--cluster 'qsub -q {cluster.qname} -l nodes={cluster.nodes}:ppn={cluster.ppn} -l mem={cluster.mem} -l walltime={cluster.time} -M {cluster.account} -m ea -o error_files/ -e error_files/' \
--use-conda \

# generate report
snakemake \--report runs/snake_ChIP_${TIME}.html
# -o {cluster.std_oe}
