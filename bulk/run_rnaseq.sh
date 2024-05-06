## Load Nextflow and Singularity environment modules
module purge
module load dangpu_libs java/11.0.15 nextflow/22.10.6 singularity/3.8.0 python/3.7.13 nf-core/2.7.2

# set up bash environment variables for memory
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_HOME=/cache/nxf-home
export NXF_TEMP=/scratch/temp/bns631
export NXF_SINGULARITY_CACHEDIR=cache/singularity-images

mkdir -p $NXF_SINGULARITY_CACHEDIR
mkdir -p $NXF_HOME
mkdir $NXF_TEMP

nextflow run nf-core/rnaseq  -r 3.10.1 -profile ku_sund_dangpu --outdir Output_skip_trimming --skip_trimming true  --input rna_samples.csv --genome hg38
