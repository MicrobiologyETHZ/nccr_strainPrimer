snakemake --use-conda -k --cluster "DIR=\$(dirname {params.qoutfile}); mkdir -p \"\${{DIR}}\"; qsub -S /bin/bash -V -cwd -o {params.qoutfile} -e {params.qerrfile} -pe smp {threads} -l h_vmem={params.mem}M" --configfile src/configs/config.yaml -p -j 10 --max-jobs-per-second 1

