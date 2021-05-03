from pathlib import Path
from Bio import SeqIO
import subprocess
import os
import pandas as pd
import re

def check_call(command):
    """
    Simple wrapper to execute check_call and catch exceptions
    :param command:
    :return:
    """
    returncode = -1
    try:
        returncode = subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        raise Exception(e)
    if returncode != 0:
        raise(Exception(f'Command: {command} failed with return code {returncode}'))
    return returncode


def write_primer3_file(probe_id, probe_seq, fileName=''):
    template = f'SEQUENCE_ID={probe_id}\n' \
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
    if fileName:
        with open(fileName, 'w') as fo:
            fo.write(template)
    return template


def find_primers_for_probe(probe):
    primer3_argument = write_primer3_file(probe.id, probe.seq)
    cmd = f'primer3_core <(printf \"{primer3_argument}\")'
    primers = subprocess.check_output(cmd, shell=True, executable='/bin/bash').decode()
    return primers


def parse_primer3(primer3_out):
    primer3_to_dict = {line.split('=')[0]: line.split('=')[1] for line in primer3_out.strip().split("\n") if len(line) > 1}
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
        final_dict = {k: 'NA' for k in keys_to_return}
        final_dict['SEQUENCE_ID'] = primer3_to_dict['SEQUENCE_ID']
        return final_dict
def find_primer_coordinates(primer3_dict, genome_file, output_dir):
    # Write fasta file
    if primer3_dict['PRIMER_LEFT_0_SEQUENCE'] == 'NA':
        return ''
    header_left = primer3_dict['SEQUENCE_ID']+"_PRIMER_LEFT"
    header_right = primer3_dict['SEQUENCE_ID']+"_PRIMER_RIGHT"
    primer_pair_fasta = Path(output_dir)/f'{primer3_dict["SEQUENCE_ID"]}.fasta'
    with open(primer_pair_fasta, 'w') as fo:
        fo.write(f'>{header_left}\n{primer3_dict["PRIMER_LEFT_0_SEQUENCE"]}\n'
                 f'>{header_right}\n{primer3_dict["PRIMER_RIGHT_0_SEQUENCE"]}\n')
    alignment = align_primers(primer_pair_fasta, genome_file)
    primer_positions = parse_sam_alignment(alignment)
    os.remove(primer_pair_fasta)
    return primer_positions


def align_primers(primer_pair_fasta, genome_file):
    # Align primer fasta file against a target file with BWA, allowing only 1 error
    # sys.stderr.write("Aligning {}\n".format(target_file))
    bwa_cmd = f"bwa index {genome_file} 2> /dev/null; " \
              f"bwa aln -n 1 {genome_file} {primer_pair_fasta} 2> /dev/null | " \
              f"bwa samse {genome_file} - {primer_pair_fasta} 2> /dev/null"
    alignment = subprocess.check_output(bwa_cmd, shell=True, executable="/bin/bash").decode()
    return alignment


class FlagException(Exception):
    pass


def parse_sam_alignment(alignment):
    alignment = alignment.strip().split("\n")
    # Parse a BWA alignment and return hits in useful format
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


def find_primers_for_strain(strain_fasta, output_dir):
    records = SeqIO.parse(strain_fasta, 'fasta')
    probe_primers = []
    for record in records:
        primer3_out = find_primers_for_probe(record)
        primers = parse_primer3(primer3_out)
        genome_file = Path(output_dir)/'genomes.fasta'
        primer_positions = find_primer_coordinates(primers, genome_file, output_dir)
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
    df = pd.DataFrame(probe_primers).set_index('SEQUENCE_ID').T
    return df


def find_all_primers(output_dir, samples):
    for sample in samples:
        strain_fasta = Path(output_dir)/f"{sample}_probes.unique.fasta"
        sampleDf = find_primers_for_strain(strain_fasta, output_dir)
        sampleDf.to_csv(Path(output_dir)/f"{sample}_final_primer_list.csv")


if __name__ == "__main__":
    samples = ['Z5841', 'Z5871', 'Z5951', 'Z5971', 'Z5991']
    probe_dir = "data/probes/"
    primer_dir = "data/primers"
    find_all_primers(probe_dir, primer_dir, samples)