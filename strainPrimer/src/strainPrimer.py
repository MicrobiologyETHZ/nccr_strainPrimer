"""
strainPrimer

Step 1: CATCH
Step 2: BLAST
Step 3: PRIMER3
Step 4: in silico PCR

"""
import argparse
import pandas as pd
import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from pathlib import Path
import os
import re
import subprocess

import logging


def get_logger(log_file_path):
    #  Create a custom logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    # Create handlers
    stream_handler = logging.StreamHandler(sys.stdout)
    file_handler = logging.FileHandler(log_file_path)
    stream_handler.setLevel(logging.INFO)
    file_handler.setLevel(logging.ERROR)
    # Create formatters and add it to handlers
    stream_format = logging.Formatter('%(asctime)s- %(levelname)s - %(message)s')
    file_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stream_handler.setFormatter(stream_format)
    file_handler.setFormatter(file_format)
    # Add handlers to the logger
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)
    return logger


def check_call(command, logger):
    """
    Simple wrapper to execute check_call and catch exceptions
    :param command: shell command to run
    :param logger: logger instance
    :return: return code of the command
    """
    try:
        returncode = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        logger.exception(e)
        raise Exception(e)
    if returncode != 0:
        logger.exception(f'Command: {command} failed with return code {returncode}')
        raise(Exception(f'Command: {command} failed with return code {returncode}'))
    return returncode


def run_catch_design(genomes, fasta_with_probes, logger, probe_length=150, genome_coverage=0.1):
    """
    Run CATCH to find unique parts for each of the genomes
    :param genomes: list of genomes to analyze
    :param fasta_with_probes: output file name in which probes will be written
    :param probe_length: length of the probes
    :param genome_coverage: % of genome that should be covered by the probes
    :param logger: logger instance
    :return: returncode of CATCH

    """
    genomes_string = ' '.join([str(genome) for genome in genomes])
    catch_cmd = f"design.py {genomes_string} -pl {probe_length} -c {genome_coverage} --identify -o {fasta_with_probes}"
    returncode = check_call(catch_cmd, logger)
    return returncode


def blast_probes_against_input_genomes(genomes, output_dir, fasta_with_probes, probe_blast_hits, logger):
    """
    Make sure CATCH identified unique regions by blasting against the genomes of interest
    :param genomes:
    :param output_dir:
    :param fasta_with_probes:
    :param probe_blast_hits:
    :param logger
    :return:
    """
    # todo runds locally without errors, test on cluster
    genomes_string = ' '.join([str(genome) for genome in genomes])
    genomes_fasta = Path(output_dir)/'genomes.fasta'
    blast_cmd = f'zcat {genomes_string} > {genomes_fasta}; ' \
                f'makeblastdb -in {genomes_fasta} -dbtype nucl ;'\
                f'blastn -db {genomes_fasta} -out {probe_blast_hits} ' \
                f'-query {fasta_with_probes} -outfmt ' \
                f'"6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand" ' \
                f'-num_threads 16'
    returncode = check_call(blast_cmd, logger)
    return returncode


def check_uniquenes(probe_blast_hits):
    df = pd.read_table(probe_blast_hits, header=None)
    df.columns = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand".split()
    uniq = df.groupby('qseqid').sseqid.count().reset_index()
    uniq = uniq[uniq.sseqid == 1].qseqid.values
    df_with_uniq_probes = df[df.qseqid.isin(uniq)]
    df_with_uniq_probes = df_with_uniq_probes[(df_with_uniq_probes.pident == 100)
                                              & (df_with_uniq_probes.length == 150)
                                              & (df_with_uniq_probes.evalue < 10e-4)]
    df_with_uniq_probes['sampleID'] = df_with_uniq_probes['sseqid'].str.split("_", expand=True)[0]
    return df_with_uniq_probes


def get_record(seq, name):
    return SeqRecord(Seq(seq), id=name, name="", description="",)


def split_by_sample(df_with_uniq_probes, output_dir="../data/blast"):
    """
    :param df_with_uniq_probes:
    :param output_dir:
    :return:
    """
    for sample, g in df_with_uniq_probes.groupby('sampleID'):
        out_file = f"{output_dir}/{sample}_probes.fasta"
        seqs = g.qseq.values
        ids = g.qseqid.values
        records = [get_record(seq, name) for seq, name in zip(seqs, ids)]
        SeqIO.write(records, out_file, "fasta")
    return list(df_with_uniq_probes.sampleID.unique())


def get_unique_probes(genomes, output_dir, force_catch, logger):
    fasta_with_probes = Path(output_dir)/"probes.fasta"
    probe_blast_hits = Path(output_dir)/"probes_blast_internal_db.txt"
    if not force_catch and fasta_with_probes.is_file():
        logger.info("Skipping CATCH")
    else:
        logger.info("Running CATCH to design unique probes") # todo add check if this step completed sucessfully
        run_catch_design(genomes, fasta_with_probes, logger)
    logger.info("Running BLAST to make sure probes are unique")
    blast_probes_against_input_genomes(genomes, output_dir, fasta_with_probes, probe_blast_hits, logger)
    logger.info("Writing unique probes for each genome")
    df_with_uniq_probes = check_uniquenes(probe_blast_hits)
    samples = split_by_sample(df_with_uniq_probes, output_dir)
    return samples

# -------------------------------- #
#    Step 2: CHECK BACKGROUND      #
# -------------------------------- #


def check_backgorund(sample_probe_fasta, ncbi_db, negative_taxid, background_blast_output, logger, max_targets=1):
    """
    BLAST each probe against NCBI db, see if there are any matches
    :param sample_probe_fasta:
    :param ncbi_db:
    :param background_blast_output:
    :param logger
    :param max_targets
    :return:
    """

    negative_taxa = f" -negative_taxids {negative_taxid} " if negative_taxid else ''
    max_target_seqs = f" -max_target_seqs {max_targets} " if max_targets > 0 else ''
    cmd = f'blastn -db {ncbi_db} -out {background_blast_output} ' \
          f'-query {sample_probe_fasta} {negative_taxa} -outfmt ' \
          f'"6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand" ' \
          f'-num_threads 32 ' \
          f'{max_target_seqs}'
    return check_call(cmd, logger)


def get_probes_with_no_background(sample_probe_fasta, background_blast_output):
    try:
        df = pd.read_table(background_blast_output, header=None)
        probes_with_background_hits = df[0].unique()
    except pd.errors.EmptyDataError:
        probes_with_background_hits = []
    finally:
        sample_probes = SeqIO.parse(sample_probe_fasta, 'fasta')
        uniq_probes = []
        for probe in sample_probes:
            if probe.id not in probes_with_background_hits:
                uniq_probes.append(probe)
        SeqIO.write(uniq_probes, Path(sample_probe_fasta).with_suffix(".unique.fasta"), "fasta")
        return Path(sample_probe_fasta).with_suffix(".unique.fasta")


def get_all_probes_without_background(samples, ncbi_db, negative_taxid, output_dir, logger):
    logger.info(f"Checking the probes against {ncbi_db}")
    unique_probe_files = []
    for sample in samples:
        logger.info(f"BLASTING {sample}")
        sample_probe_fasta = Path(output_dir)/f'{sample}_probes.fasta'
        background_blast_output = Path(output_dir)/f'{sample}.checkBackground.blast'
        check_backgorund(sample_probe_fasta, ncbi_db, negative_taxid, background_blast_output, logger)
        logger.info(f"BLASTING {sample} done, writing unique probes")
        unique_fasta = get_probes_with_no_background(sample_probe_fasta, background_blast_output)
        unique_probe_files.append(unique_fasta)
    return unique_probe_files


# -------------------------------- #
#      Step 3: Run Primer3         #
# -------------------------------- #

# todo need to parallelize this to make it faster !!!

def get_primer3_argument(probe_id, probe_seq):
    return f'SEQUENCE_ID={probe_id}\n' \
               f'SEQUENCE_TEMPLATE={probe_seq}\n' \
               f'PRIMER_TASK=generic\n' \
               f'PRIMER_PICK_LEFT_PRIMER=1\n' \
               f'PRIMER_PICK_INTERNAL_OLIGO=0\n' \
               f'PRIMER_PICK_RIGHT_PRIMER=1\n' \
               f'PRIMER_OPT_SIZE=18\n' \
               f'PRIMER_MIN_SIZE=15\n' \
               f'PRIMER_MAX_SIZE=21\n' \
               f'PRIMER_MAX_NS_ACCEPTED=1\n' \
               f'PRIMER_PRODUCT_SIZE_RANGE=120-150\n' \
               f'P3_FILE_FLAG=1\n' \
               f'PRIMER_EXPLAIN_FLAG=1\n' \
               f'=\n'


def find_primers_for_probe(probe):
    """
    :param probe: SeqIO object
    :return:
    """
    primer3_argument = get_primer3_argument(probe.id, probe.seq)
    primer3_command = f'primer3_core <(printf \"{primer3_argument}\")'
    primer3_output = subprocess.check_output(primer3_command, shell=True, executable='/bin/bash').decode()
    return primer3_output


def parse_primer3_output(primer3_output):
    primer3_to_dict = {line.split('=')[0]: line.split('=')[1]
                       for line in primer3_output.strip().split("\n") if len(line) > 1}
    keys_to_return = ['SEQUENCE_ID', 'PRIMER_PAIR_0_PENALTY', 'PRIMER_LEFT_0_PENALTY', 'PRIMER_RIGHT_0_PENALTY',
                      'PRIMER_LEFT_0_SEQUENCE', 'PRIMER_RIGHT_0_SEQUENCE', 'PRIMER_LEFT_0_TM', 'PRIMER_RIGHT_0_TM',
                      'PRIMER_LEFT_0_GC_PERCENT', 'PRIMER_RIGHT_0_GC_PERCENT', 'PRIMER_LEFT_0_SELF_ANY_TH',
                      'PRIMER_RIGHT_0_SELF_ANY_TH', 'PRIMER_LEFT_0_SELF_END_TH', 'PRIMER_RIGHT_0_SELF_END_TH',
                      'PRIMER_LEFT_0_HAIRPIN_TH', 'PRIMER_RIGHT_0_HAIRPIN_TH', 'PRIMER_LEFT_0_END_STABILITY',
                      'PRIMER_RIGHT_0_END_STABILITY', 'PRIMER_PAIR_0_COMPL_ANY_TH', 'PRIMER_PAIR_0_COMPL_END_TH',
                      'PRIMER_PAIR_0_PRODUCT_SIZE']
    if int(primer3_to_dict['PRIMER_PAIR_NUM_RETURNED']) > 0:
        return {k: primer3_to_dict[k] for k in keys_to_return}
    else:
        primer_dict = {k: 'NA' for k in keys_to_return}
        primer_dict['SEQUENCE_ID'] = primer3_to_dict['SEQUENCE_ID']
        return primer_dict


def align_primers(primer_pair_fasta, genome_file, logger):
    """
    Align primer fasta file against a target file with BWA, allowing only 1 error
    :param primer_pair_fasta:
    :param genome_file:
    :param logger
    :return:
    """
    logger.info(f'Aligning {Path(primer_pair_fasta).name} to {Path(genome_file).name}')
    bwa_cmd = f"bwa index {genome_file} 2> /dev/null; " \
              f"bwa aln -n 1 {genome_file} {primer_pair_fasta} 2> /dev/null | " \
              f"bwa samse {genome_file} - {primer_pair_fasta} 2> /dev/null"
    alignment = subprocess.check_output(bwa_cmd, shell=True, executable="/bin/bash").decode()
    return alignment


class FlagException(Exception):
    pass


def parse_sam_alignment(alignment):
    """
    Parse a BWA alignment and return hits in useful format
    :param alignment:
    :return:
    """
    alignment = alignment.strip().split("\n")
    header_marker = '@'
    unmapped_flag = '4'
    mapped_forward = '0'
    mapped_reverse = '16'
    primer_positions = {}
    for line in alignment:
        if not line.startswith(header_marker):
            fields = line.strip().split("\t")
            primer_id = fields[0]
            contig_id = fields[2]
            start_position = fields[3]
            primer_sequence = fields[9]
            if fields[1] == unmapped_flag:
                continue
            elif fields[1] == mapped_forward:
                strand = '+'
            elif fields[1] == mapped_reverse:
                strand = '-'
            else:
                raise FlagException(f"Cannot interprete FLAG: {fields[1]}")
            edit_distance = re.sub("NM:i:", "", fields[12])
            if primer_id not in primer_positions.keys():  # todo refactor with get
                primer_positions[primer_id] = [(primer_sequence, contig_id, strand, start_position, edit_distance)]
            else:
                primer_positions[primer_id] += [(primer_sequence, contig_id, strand, start_position, edit_distance)]
    return primer_positions


def find_primer_coordinates(primer_dict, genome_file, output_dir, logger):
    """

    :param primer_dict:
    :param genome_file:
    :param output_dir:
    :param logger
    :return:
    """
    if primer_dict['PRIMER_LEFT_0_SEQUENCE'] == 'NA':
        return ''
    header_left = primer_dict['SEQUENCE_ID'] + "_PRIMER_LEFT"
    header_right = primer_dict['SEQUENCE_ID'] + "_PRIMER_RIGHT"
    primer_pair_fasta = Path(output_dir)/f'{primer_dict["SEQUENCE_ID"]}.fasta'
    with open(primer_pair_fasta, 'w') as fo:
        fo.write(f'>{header_left}\n{primer_dict["PRIMER_LEFT_0_SEQUENCE"]}\n'
                 f'>{header_right}\n{primer_dict["PRIMER_RIGHT_0_SEQUENCE"]}\n')
    alignment = align_primers(primer_pair_fasta, genome_file, logger)
    primer_positions = parse_sam_alignment(alignment)
    os.remove(primer_pair_fasta)
    return primer_positions


def find_primers_for_strain(strain_fasta, output_dir, logger):
    records = SeqIO.parse(strain_fasta, 'fasta')
    probe_primers = []
    for record in records:
        logger.info(f"Finding primers for {record.id}")
        logger.info('Running Primer3')
        primer3_out = find_primers_for_probe(record)
        primers = parse_primer3_output(primer3_out)
        genome_file = Path(output_dir)/'genomes.fasta'  # todo don't hard code this?
        # todo align to each genome separately?
        logger.info("Finding primer positions in the genome")
        primer_positions = find_primer_coordinates(primers, genome_file, output_dir, logger)
        if primer_positions:
            left_key = primers['SEQUENCE_ID']+"_PRIMER_LEFT"
            right_key = primers['SEQUENCE_ID'] + "_PRIMER_RIGHT"
            primers['PRIMER_LEFT_0_SEQID'] = [position[1] for position in primer_positions[left_key]]
            primers['PRIMER_LEFT_0_POSITION'] = [position[3] for position in primer_positions[left_key]]
            primers['PRIMER_LEFT_0_STRAND'] = [position[2] for position in primer_positions[left_key]]
            primers['PRIMER_RIGHT_0_SEQID'] = [position[1] for position in primer_positions[right_key]]
            primers['PRIMER_RIGHT_0_POSITION'] = [position[3] for position in primer_positions[right_key]]
            primers['PRIMER_RIGHT_0_STRAND'] = [position[2] for position in primer_positions[right_key]]
        probe_primers.append(primers)
        logger.info("Done")
    logger.info(f"Finished {Path(strain_fasta).stem}")
    if probe_primers:
        df = pd.DataFrame(probe_primers).set_index('SEQUENCE_ID').T
        return df
    else:
        return pd.DataFrame()


def find_all_primers(output_dir, samples, logger):
    for sample in samples:
        logger.info(f"Finding pimers for {sample}")
        strain_fasta = Path(output_dir)/f"{sample}_probes.unique.fasta"
        sampleDf = find_primers_for_strain(strain_fasta, output_dir, logger)
        if sampleDf.empty:
            logger.warning("No unique primers found")
        else:
            logger.info(f"Writing results for {sample} to file")
            sampleDf.to_csv(Path(output_dir)/f"{sample}_final_primer_list.csv")


def get_final_probes_fasta(sample, output_dir, logger):
    output_dir = Path(output_dir)
    logger.info(f"Checking uniqueness of probes for {sample}")
    if Path(output_dir/f"{sample}_final_primer_list.csv").is_file() and Path(output_dir/f"{sample}_probes.fasta").is_file():
        probes = pd.read_csv(Path(output_dir/f"{sample}_final_primer_list.csv")).dropna(axis=1).columns
        out_file = f"{output_dir}/{sample}_final_probes.fasta"
        records = SeqIO.parse(Path(output_dir)/f"{sample}_probes.fasta", 'fasta')
        final_records = [record for record in records if record.id in probes]
        SeqIO.write(final_records, out_file, "fasta")
        return out_file
    elif not Path(output_dir/f"{sample}_probes.fasta").is_file():
        logger.info(f"No probes found for {sample}")
        return None
    else:
        logger.info(f"No final primer list found for {sample}")
        return None


def get_background_genomes_with_final_probes(sample, final_probes_fasta, ncbi_db, output_dir, logger):
    # blast them against whole NCBI
    if final_probes_fasta:
        background_blast_output = Path(output_dir)/ f'{sample}.finalProbes.blast'
        check_backgorund(sample_probe_fasta=final_probes_fasta, ncbi_db=ncbi_db,
                         negative_taxid=None, background_blast_output=background_blast_output,
                         logger=logger, max_targets=0)

    else:
        logger.info("No final probes fasta")
        return None




# Run the whole pipeline

def get_primers(genome_dir, output_dir, ncbi_db, negative_taxid, force_catch):
    genomes = [str(genome) for genome in Path(genome_dir).iterdir()]
    logger = get_logger(Path(output_dir)/"strainPrimer.log")
    logger.info('Step 1')
    # Step 1: Identifies unique probes and returns a list of sample names (for each genome analysed)
    samples = get_unique_probes(genomes, output_dir, force_catch, logger)
    logger.info(f'Finished Step 1. These are the samples analysed: {samples}')
    logger.warning('Step 2')
    # Step 2: BLAST
    get_all_probes_without_background(samples, ncbi_db, negative_taxid, output_dir, logger)
    # Step 3: Primer3
    logger.info('Step 3')
    find_all_primers(output_dir, samples, logger)
    logger.info("All Done!")
    return None

def parse_args():
    parser = argparse.ArgumentParser(
        description='Find unique probes in genomes and design qPCR primers for them')
    parser.add_argument('-g', '--genomeDir',   help='Directory containing genomes in fasta format')
    parser.add_argument('-f', '--force_catch', action='store_true', help='Force CATCH to run, even if probes.fasta already exists')
    parser.add_argument('-db', '--database',  default='/nfs/cds/Databases/BLAST-NCBI/nt/June_2019/nt_v5',
                        help='BLAST reference database')
    parser.add_argument('-o', '--outDir',  help='Output directory')
    parser.add_argument('-e,', '--exclude', help='Taxid to exclude') # todo make sure can handle multiple

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    get_primers(args.genomeDir, args.outDir, args.database, args.exclude, args.force_catch)


