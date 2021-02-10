"""

Looks for unique blast hits for the probes and sorts them by sample they map to
Writes out fasta file for each sample

"""
import pandas as pd
import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

def get_unique(blast_file):
    df = pd.read_table(blast_file, header=None)
    df.columns = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sstrand".split()
    uniq = df.groupby('qseqid').sseqid.count().reset_index()
    uniq = uniq[uniq.sseqid == 1].qseqid.values
    df = df[df.qseqid.isin(uniq)]
    df = df[(df.pident==100) & (df.length==150) & (df.evalue<10e-4)]
    df['sampleID'] = df['sseqid'].str.split("_", expand=True)[0]
    return df


def split_by_sample(df, out_dir="../data/blast"):
    for sample, g in df.groupby('sampleID'):
        out_file = f"{out_dir}/{sample}_probes.fasta"
        seqs = g.qseq.values
        ids = g.qseqid.values
        records = [get_record(seq, name) for seq, name in zip(seqs, ids)]
        SeqIO.write(records, out_file, "fasta")


def get_record(seq, name):
    return SeqRecord(Seq(seq), id=name, name="",description="",)


def main():
    blast_file = sys.argv[1]
    out_dir = sys.argv[2]
    df = get_unique(blast_file)
    split_by_sample(df, out_dir)


if __name__ == "__main__":
    main()