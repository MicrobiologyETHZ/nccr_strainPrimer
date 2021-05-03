from pathlib import Path
from strainPrimer.src.get_unique_probes import *
GENOMEDIR = "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/mkmherzog/hfd/scratch/genomes_for_strainPrimer_test"
OUTDIR = "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/mkmherzog/hfd/scratch/27_04_strainPrimer"

def test_run_catch_design():
    genomes = [str(genome) for genome in Path(GENOMEDIR).iterdir()]
    print(" ".join(genomes))
    fasta_with_probes = Path(OUTDIR)/"probes.fasta"
    run_catch_design(genomes, fasta_with_probes)


#test_run_catch_design()

def test_blast_probes_against_input_genomes():
    genomes = [str(genome) for genome in Path(GENOMEDIR).iterdir()]
    fasta_with_probes = Path(OUTDIR)/"probes.fasta"
    probe_blast_hits = Path(OUTDIR)/"probes_blast_internal_db.txt"
    blast_probes_against_input_genomes(genomes, OUTDIR, fasta_with_probes, probe_blast_hits)

#test_blast_probes_against_input_genomes()

def test_get_unique_probes():
    genomes = [str(genome) for genome in Path(GENOMEDIR).iterdir()]
    samples = get_unique_probes(genomes, OUTDIR)
    print(samples)

#test_get_unique_probes()
def test_check_background():
    sample_probe_fasta = Path(OUTDIR)/'Z6011_probes.fasta'
    background_blast_output = Path(OUTDIR)/'Z6011.checkBackground.blast'
    NCBI_db = '/nfs/cds/Databases/BLAST-NCBI/nt/June_2019/nt_v5'
    check_backgorund(sample_probe_fasta, NCBI_db, background_blast_output)
    unique_hits(sample_probe_fasta, background_blast_output)

test_check_background()