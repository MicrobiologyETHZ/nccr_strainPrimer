import click
from pathlib import Path
from strainPrimer.src.strainPrimer import get_primers
@click.group()
def main():
    pass

# Get unique probes
@main.command(help='Run CATCH and BLAST to get unique probes ')
@click.option('--output_dir', '-o',  help='Output Directory')
@click.option('--ncbi_db', '-db',  help='NCBI BLAST db')
@click.option('--negative_taxid', '-ntax',  help='taxid to mask, ex. 562 for E.coli')
@click.option('--genome_dir', '-tn', help='genome directory')
def probes(genome_dir, ncbi_db, negative_taxid, output_dir):
    get_primers(genome_dir, output_dir, ncbi_db, negative_taxid, )
