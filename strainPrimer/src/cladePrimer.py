import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

def parse_cdhit_output(cdhit_clstr_file):
    clusters = []
    with open(cdhit_clstr_file, 'r') as fh:
        cluster = ''
        for line in fh.readlines():
            c = []
            if line.startswith(">"):
                cluster = line.strip('>\n')
            else:
                members = line.split('\t')
                member_ix = members[0]
                member_length = int(members[1].split(',')[0].strip('a'))
                member_id = members[1].split(',')[1].strip(' >').split("...")[0]
                member_ident = members[1].split('...')[1].strip('at% \n')
                member_ident = 100 if member_ident == '*' else float(member_ident)
                c = [cluster, member_ix, member_length, member_id, member_ident]
                clusters.append(c)
    df = pd.DataFrame(clusters, columns=['cluster', 'cluster_ix', 'prot_length', 'gene_id', 'percent_ident'])
    return df


def get_cluster_nucl_sequence(concat_fasta, cluster_genes, out_file):
    gene_records = []
    records = SeqIO.parse(concat_fasta, 'fasta')
    for record in records:
        if record.id in cluster_genes:
            gene_records.append(record)
    SeqIO.write(gene_records, out_file, "fasta")


def alilgnment_to_df(alignment_file):
    rids = []
    seqs = []
    records = SeqIO.parse(alignment_file, 'fasta')
    for record in records:
        rids.append(record.id)
        seqs.append(str(record.seq))
    alignments = [rids, seqs]

    df = pd.DataFrame(alignments, index =['gene_id', 'seq']).T
    positions = pd.DataFrame(df.seq.apply(list).to_list())
    return pd.concat([df['gene_id'], positions], axis=1)


if __name__ == "__main__":
    df = alilgnment_to_df("/Users/ansintsova/git_repos/strainPrimer/strainPrimer/tests/test_files/c5930.strain.alignment")
    print(df.head())
    # cdhit_clstr_file = "/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/mkmherzog/hfd/scratch/assembly/protein_hfd.clstr"
    # df = parse_cdhit_output(cdhit_clstr_file)
    # df.to_csv("/nfs/nas22/fs2202/biol_micro_bioinf_nccr/hardt/mkmherzog/hfd/scratch/assembly/protein_hfd.csv")
 #    cluster_genes = ['Z5741_00150',
 # 'Z5761_00356',
 # 'Z5761_04232',
 # 'Z5781_04834',
 # 'Z5801_01688',
 # 'Z5821_01594',
 # 'Z5821_02825',
 # 'Z5841_00442',
 # 'Z5871_03227',
 # 'Z5894_00540',
 # 'Z5894_04209',
 # 'Z5931_00130',
 # 'Z5931_04132',
 # 'Z5951_01069',
 # 'Z5971_01645',
 # 'Z5991_04494',
 # 'Z6011_00542',
 # 'Z6011_03865',
 # 'Z6021_03660',
 # 'Z6021_03684',
 # 'Z6023_00612',
 # 'Z6024_01161',
 # 'Z6024_04551',
 # 'Z6025_02578',
 # 'Z6026_01302',
 # 'Z6027_01862',
 # 'Z6027_02754']
 #    get_cluster_nucl_sequence("/Users/ansintsova/git_repos/strainPrimer/strainPrimer/tests/test_files/hfd.strain.ffn",
 #                              cluster_genes, "/Users/ansintsova/git_repos/strainPrimer/strainPrimer/tests/test_files/c5930.strain.ffn")