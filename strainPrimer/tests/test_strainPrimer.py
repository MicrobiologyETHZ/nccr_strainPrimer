from pathlib import Path
import logging
from strainPrimer.src.strainPrimer import blast_probes_against_input_genomes


def test_blast_probes_against_input_genomes():
    root = Path("/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/mkmherzog/hfd/scratch/genomes_for_strainPrimer_test")
    genomes = [root/"Z6011.scaffolds.min500.fasta.gz",
               root/"Z6021.scaffolds.min500.fasta.gz",
               root/"Z6023.scaffolds.min500.fasta.gz"]
    output_dir = root
    fasta_with_probes = root/"test_probes.fasta"
    probe_blast_hits = root/'test_blast_hits.out'
    logger = logging.Logger(name='Test')
    rc = blast_probes_against_input_genomes(genomes, output_dir, fasta_with_probes, probe_blast_hits, logger)
    print(rc)

test_blast_probes_against_input_genomes()
