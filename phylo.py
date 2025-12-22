#O código não é reprodutível na sua forma atual, pois foi inicialmente desenvolvido e testado no Google Colab, onde é possível executar comandos de sistema diretamente no ambiente de trabalho. Nesse contexto, foi possível ligar o alinhamento múltiplo de sequências (MSA) à construção da árvore filogenética num único fluxo de execução.
#Ao tentar executar o mesmo código localmente, surgiram dificuldades na instalação e integração de ferramentas externas, nomeadamente o MAFFT, o que impediu a execução completa do pipeline dentro de um único script. Como resultado, o alinhamento teve de ser realizado externamente, sendo executada localmente apenas uma versão simplificada do código, suficiente para gerar a árvore filogenética a partir de um alinhamento já existente.
#A resolução destes problemas encontra-se ainda em desenvolvimento, com o objetivo de tornar o pipeline totalmente reprodutível num ambiente local.


from pathlib import Path
import shutil
import subprocess

INPUT_FASTA = "blast_hits.fasta"  
USE_GDOWN = False
URL_FASTA = "https://drive.google.com/uc?id=1bl6HxzzBHjfObybojQOpPjr51mizVqjc"

MAFFT_EXE = shutil.which("mafft") or "mafft"

MSA_FASTA = "hits.msa.fasta"
PHYLIP_OUT = "msa.phy"
LEGEND_OUT = "legend.txt"
TREE_UPGMA_NWK = "tree_upgma.nwk"
TREE_NJ_NWK = "tree_nj.nwk"


from Bio import SeqIO

seqs = list(SeqIO.parse(INPUT_FASTA, "fasta"))
assert len(seqs) >= 3
assert all(len(s.seq) > 0 for s in seqs)

print("Nº seqs:", len(seqs))
print("Primeiro ID:", seqs[0].id)
print("Tamanho 1ª seq:", len(seqs[0].seq))

from pathlib import Path

INPUT_FASTA = "blast_hits.fasta"  
MSA_FASTA = "hits.msa.fasta"
TREE_UPGMA_NWK = "tree_upgma.nwk"
TREE_NJ_NWK = "tree_nj.nwk"

assert Path(INPUT_FASTA).exists(), f"Não encontrei {INPUT_FASTA} na pasta do notebook"
print("INPUT_FASTA:", INPUT_FASTA)

from Bio import SeqIO
seqs = list(SeqIO.parse(INPUT_FASTA, "fasta"))
assert len(seqs) >= 3

from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.5

def merge_into_msa(msa_seqs, ref_new, seq_new):
    old_ref = msa_seqs[0]

    i_old = 0
    i_new = 0
    out_msa = [""] * len(msa_seqs)
    out_seq = ""

    while i_old < len(old_ref) or i_new < len(ref_new):
        c_old = old_ref[i_old] if i_old < len(old_ref) else None
        c_new = ref_new[i_new] if i_new < len(ref_new) else None

        if c_old == "-" and c_new == "-":
            for k in range(len(msa_seqs)):
                out_msa[k] += "-"
            out_seq += "-"
            i_old += 1
            i_new += 1

        elif c_old == "-":
            for k in range(len(msa_seqs)):
                out_msa[k] += "-"
            out_seq += "-"   
            i_old += 1

        elif c_new == "-":
            for k in range(len(msa_seqs)):
                out_msa[k] += msa_seqs[k][i_old]
            out_seq += "-"
            i_new += 1
            i_old += 1

        else:
            for k in range(len(msa_seqs)):
                out_msa[k] += msa_seqs[k][i_old]
            out_seq += seq_new[i_new]
            i_old += 1
            i_new += 1

    return out_msa, out_seq

msa_ids = [seqs[0].id]
msa_aligned = [str(seqs[0].seq)]

for rec in seqs[1:]:
    ref = msa_aligned[0].replace("-", "")
    aln = aligner.align(ref, str(rec.seq))[0]
    ref_aln, seq_aln = str(aln[0]), str(aln[1])
    msa_aligned, new_aligned = merge_into_msa(msa_aligned, ref_aln, seq_aln)
    msa_ids.append(rec.id)
    msa_aligned.append(new_aligned)

msa_records = [
    SeqRecord(Seq(msa_aligned[i]), id=msa_ids[i], description="")
    for i in range(len(msa_ids))
]
SeqIO.write(msa_records, MSA_FASTA, "fasta")
print("MSA criado (sem executáveis):", MSA_FASTA)

from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

alignment = AlignIO.read(MSA_FASTA, "fasta")

dm = DistanceCalculator("identity").get_distance(alignment)
constructor = DistanceTreeConstructor()

tree_upgma = constructor.upgma(dm)
tree_nj = constructor.nj(dm)

Phylo.write(tree_upgma, TREE_UPGMA_NWK, "newick")
Phylo.write(tree_nj, TREE_NJ_NWK, "newick")

print("UPGMA:")
Phylo.draw_ascii(tree_upgma)
print("\nNJ:")
Phylo.draw_ascii(tree_nj)
