from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import subprocess
import os
import matplotlib.pyplot as plt
import re


def msa_from_blastp_xml(
    xml_file: str,
    query_fasta: str,
    output_fasta: str = "hits_prot.fasta",
    n_hits: int = 10,
    email: str = "teu_email@email.com",
):
    """
    Os resultados BLASTp em formato XML foram processados, sendo selecionados os n melhores alinhamentos.
    Para cada hit identificado, foi obtida a sequência proteica completa associada ao respetivo accession,
    recorrendo à base de dados protein do NCBI. A sequência proteica da query foi incluída juntamente com
    as proteínas homólogas recuperadas, formando um ficheiro FASTA multi-sequência. Este ficheiro serviu
    de base para o alinhamento múltiplo de sequências e para a inferência filogenética subsequente,
    baseada em sequências de aminoácidos.
    """

    Entrez.email = email

    records = []
    id_to_organism = {}

    print("Proteínas incluídas:")

    with open(xml_file) as f:
        blast_record = NCBIXML.read(f)

    match = re.search(r"blast_result_(\d+)_blastp", xml_file)
    query_number = match.group(1) if match else "unknown"

    query = SeqIO.read(query_fasta, "fasta")

    query_id = "query"
    query.id = query_id
    query.name = query_id
    query.description = f"Query_{query_number}"

    records.append(query)
    id_to_organism[query_id] = query.description
    print(f" - {query.description}")



    if not blast_record.alignments:
        raise ValueError("O XML não contém alinhamentos BLASTp.")

    for aln in blast_record.alignments[:n_hits]:

        organism = aln.hit_def
        organism = organism.split(">")[0].strip()

        handle = Entrez.efetch(
            db="protein",
            id=aln.accession,
            rettype="fasta",
            retmode="text"
        )

        rec = SeqIO.read(handle, "fasta")
        handle.close()

        safe_id = f"prot{len(records)}"
        rec.id = safe_id
        rec.name = safe_id
        rec.description = organism

        records.append(rec)
        id_to_organism[safe_id] = organism

        print(f" - {organism}")

    SeqIO.write(records, output_fasta, "fasta")
    print(f"\nFASTA criado: {output_fasta}")
    print(f"Total de proteínas: {len(records)}")

    return output_fasta, id_to_organism


def run_mafft(
    input_fasta: str,
    mafft_path: str,
    output_fasta: str = "hits_prot_aligned.fasta"
):
    """
    O alinhamento múltiplo de sequências proteicas foi realizado utilizando o software MAFFT.
    O ficheiro FASTA multi-sequência, contendo a proteína query e as proteínas homólogas
    recuperadas por BLASTp, foi utilizado como input, sendo aplicado o modo --auto, que
    permite ao MAFFT selecionar automaticamente a estratégia de alinhamento mais adequada.
    O alinhamento resultante foi guardado num novo ficheiro FASTA, que serviu de base para
    a análise filogenética subsequente.
    """

    print("\nA correr MAFFT...")

    with open(output_fasta, "w") as out:
        subprocess.run(
            [mafft_path, "--auto", input_fasta],
            stdout=out,
            stderr=subprocess.DEVNULL,
            check=True,
            shell=True
        )

    print(f"MSA concluído: {output_fasta}")
    return output_fasta


def phylo_tree_prot(aligned_fasta: str, id_to_organism: dict):
    """
    A função phylo_tree_prot constrói uma árvore filogenética a partir de um alinhamento
    múltiplo de sequências proteicas. O alinhamento é utilizado para calcular uma matriz
    de distâncias baseada na identidade entre aminoácidos, refletindo a proporção de
    diferenças observadas entre as sequências. Com base nesta matriz de distâncias,
    a árvore filogenética é inferida recorrendo ao método UPGMA e posteriormente
    visualizada de forma gráfica, permitindo a interpretação das relações evolutivas
    entre as proteínas analisadas.
    """

    alignment = AlignIO.read(aligned_fasta, "fasta")

    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    for clade in tree.get_nonterminals():
        clade.name = None

    fig, ax = plt.subplots(figsize=(20, 10))  
    plt.rcParams.update({"font.size": 14})

    Phylo.draw(
        tree,
        axes=ax,
        label_func=lambda clade: id_to_organism.get(clade.name)
        if clade.is_terminal()
        else None,
        do_show=False
    )

    x_left, x_right = ax.get_xlim()
    ax.set_xlim(x_left, x_right * 2.0)

    plt.subplots_adjust(left=0.05, right=0.72)

    plt.show()


    return tree


if __name__ == "__main__":

    xml_file = "...blastp.xml"
    query_fasta = "...gene.fasta"
    email = "...@email.com"

    mafft_path = r"C:\Users\..."
    # O caminho para o executável do MAFFT deve ser adaptado ao sistema operativo e à instalação local do utilizador que executar o código.

    fasta_hits, id_to_organism = msa_from_blastp_xml(
        xml_file=xml_file,
        query_fasta=query_fasta,
        email=email,
        n_hits=10
    )

    aligned = run_mafft(fasta_hits, mafft_path)

    tree = phylo_tree_prot(aligned, id_to_organism)
