module load nextflow/23.10.0
bsub -W 56:00 -R "rusage[mem=16G]" "nextflow run main.nf -c $1 -resume"
