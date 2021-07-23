from pathlib import Path
import os
import re
import subprocess


def align_primers(primer_pair_fasta, genome_file):
    """
    Align primer fasta file against a target file with BWA, allowing only 1 error
    :param primer_pair_fasta:
    :param genome_file:
    :param logger
    :return:
    """

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



def find_primer_coordinates(left, right, genome_file, output_dir):
    """
    :param genome_file:
    :param output_dir:

    :return:
    """

    header_left = "_PRIMER_LEFT"
    header_right = "_PRIMER_RIGHT"
    primer_pair_fasta = Path(output_dir)/f'test.fasta'
    with open(primer_pair_fasta, 'w') as fo:
        fo.write(f'>{header_left}\n{left}\n'
                 f'>{header_right}\n{right}\n')
    alignment = align_primers(primer_pair_fasta, genome_file)
    primer_positions = parse_sam_alignment(alignment)
    os.remove(primer_pair_fasta)
    return primer_positions


if __name__ == '__main__':
    left = 'AGAGTAACGGAGGAGCACGA'
    right = 'CAACCCGAACACCAGTGATG'
    genomes = [c for c in Path('/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/mkmherzog/hfd/scratch/genomes_for_camplicon').iterdir() if c.name.endswith('fasta')]
    print(genomes)
    output_dir = '/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/mkmherzog/hfd/scratch/30_04_camplicon'
    for genome in genomes:
        print(genome)
        print(find_primer_coordinates(left, right, genome, output_dir))