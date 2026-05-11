#PBS -S /bin/bash
#PBS -N mageckpathway
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=8:mem=24gb
#PBS -o /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/oltz_CRISPR/OU
#PBS -e /rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/oltz_CRISPR/ER

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

DIR=/rds/general/project/lms-spivakov-analysis/live/HRJ_monocytes/hILCs/oltz_CRISPR

# Run for genes identified by multiCOGS in allTraits 
mageck pathway --gene-ranking ${DIR}/biorxiv_Brown2025_suppl/geneLevel_analyses_biorxiv_screen_uppercase_COGS_ILC3_genes.txt \
        --gmt-file ${DIR}/allTraits_COGS.gmt \
        --method gsea \
        --ranking-column 2 \
        --ranking-column-2 9 \
        --keep-tmp \
        --permutation 100000 \
        -n screen3_GSEA_multiCOGS_restricted_allTraits


# the negative score is in column 3 (called by "2")
# the positive score is in column 10 (called by "9")


conda deactivate
