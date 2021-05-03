from pathlib import Path
import logging
from strainPrimer.src.strainPrimer import *
GENOMEDIR = "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/mkmherzog/hfd/scratch/genomes_for_strainPrimer_test"
OUTDIR = "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/mkmherzog/hfd/scratch/27_04_strainPrimer"


def test_find_primers_for_probe():
    strain_fasta = Path(OUTDIR) / "Z6011_probes.test.fasta"
    records = SeqIO.parse(strain_fasta, 'fasta')
    probe = list(records)[0]
    print(find_primers_for_probe(probe))

#test_find_primers_for_probe()

# def test_parse_primer3():
#     strain_fasta = Path(OUTDIR) / "Z6011_probes.test.fasta"
#     records = SeqIO.parse(strain_fasta, 'fasta')
#     probe = list(records)[0]
#     primer3_out = find_primers_for_probe(probe)
#     print(parse_primer3(primer3_out))

#test_parse_primer3()


def test_find_primers_for_strain():
    strain_fasta = Path(OUTDIR)/"Z6011_probes.test2.fasta"
    primers = find_primers_for_strain(strain_fasta, OUTDIR)
    #print(primers[0])
    print(primers.head())

#test_find_primers_for_strain()



def test_find_primers_for_strain():
    logger = logging.getLogger()
    strain_fasta = Path(OUTDIR)/'Z6021_probes.unique.fasta'
    find_primers_for_strain(strain_fasta, OUTDIR, logger)

test_find_primers_for_strain()


# def test_find_primer_coordinates():
#     strain_fasta = Path(OUTDIR) / "Z6011_probes.test.fasta"
#     file_prefix = Path(OUTDIR) / "Z6011_test"
#     primers = find_primers_for_strain(strain_fasta, file_prefix)
#     primer3_dict = primers[0]
#     genome_file = Path(OUTDIR)/"genomes.fasta"
#     alignment = find_primer_coordinates(primer3_dict, genome_file, OUTDIR)
#     print(parse_sam_alignment(alignment))
# test_find_primer_coordinates()