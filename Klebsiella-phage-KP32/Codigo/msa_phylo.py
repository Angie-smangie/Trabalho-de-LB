from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import subprocess
import os
import matplotlib.pyplot as plt
import re

def msa_from_blastn_xml(
    xml_file: str,
    query_fasta: str,
    output_fasta: str = "hits_nucl.fasta",
    n_hits: int = 10,
    email: str = "teu_email@email.com",
):
    """
    Os resultados BLASTn em formato XML foram processados, sendo  selecionados os n melhores alinhamentos. 
    Para cada hit, foi considerado apenas o HSP com maior score, correspondente à região de maior similaridade 
    entre a sequência query e a sequência da base de dados. Com base nas coordenadas desse alinhamento, foi
    extraída exclusivamente a região homóloga de cada sequência alvo, recorrendo à base de dados nuccore do NCBI.
    As regiões obtidas foram reunidas num ficheiro FASTA multi-sequência, que serviu de base para o alinhamento 
    múltiplo e para a análise filogenética subsequente.
    """
                                                                    
    Entrez.email = email

    records = []
    id_to_organism = {}

    print("Sequências incluídas:")

    query = SeqIO.read(query_fasta, "fasta")

    with open(xml_file) as f:
        blast_record = NCBIXML.read(f)

    query_label = blast_record.query

    if not query_label or query_label.lower().startswith("no definition"):
        xml_base = os.path.splitext(os.path.basename(xml_file))[0]
        gene_id = xml_base.replace("blast_result_", "")
        query_label = f"Query_GeneID_{gene_id}"

    query_id = "query"
    query.id = query_id
    query.name = query_id
    query.description = query_label

    records.append(query)
    id_to_organism[query_id] = query_label

    print(f" - {query_label}")


    if not blast_record.alignments:
        raise ValueError("O XML não tem alinhamentos BLAST.")

    for aln in blast_record.alignments[:n_hits]:
        if not aln.hsps:
            continue

        organism = aln.hit_def
        organism = re.sub(r"^MAG\s+[^:]+:\s*", "", organism)
        organism = re.sub(r",\s*(complete|partial)\s+genome", "", organism)
        organism = organism.split(">")[0].strip()

        hsp = aln.hsps[0]
        start = min(hsp.sbjct_start, hsp.sbjct_end)
        end = max(hsp.sbjct_start, hsp.sbjct_end)
        strand = 1 if hsp.sbjct_start <= hsp.sbjct_end else 2

        handle = Entrez.efetch(
            db="nuccore",
            id=aln.accession,
            rettype="fasta",
            retmode="text",
            seq_start=start,
            seq_stop=end,
            strand=strand,
        )

        rec = SeqIO.read(handle, "fasta")
        handle.close()

        safe_id = f"org{len(records)}"

        rec.id = safe_id
        rec.name = safe_id
        rec.description = organism

        records.append(rec)
        id_to_organism[safe_id] = organism

        print(f" - {organism}")

    SeqIO.write(records, output_fasta, "fasta")
    print(f"\nFASTA criado: {output_fasta}")
    print(f"Total de sequências: {len(records)}")

    return output_fasta, id_to_organism

def run_mafft(input_fasta: str, mafft_path: str, output_fasta: str = "hits_nucl_aligned.fasta"):
    
    """
    O alinhamento múltiplo de sequências foi realizado recorrendo ao software MAFFT. O ficheiro FASTA 
    multi-sequência, contendo as regiões nucleotídicas homólogas previamente extraídas, foi utilizado 
    como input, sendo aplicado o modo --auto, que permite ao MAFFT selecionar automaticamente a estratégia 
    de alinhamento mais adequada. O alinhamento resultante foi guardado num novo ficheiro FASTA, que 
    serviu de base para a análise filogenética subsequente.
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

def phylo_tree_nucl(aligned_fasta: str, id_to_organism: dict):

    """
    A função phylo_tree_nucl constrói uma árvore filogenética a partir de um alinhamento múltiplo de sequências
    nucleotídicas. O alinhamento é utilizado para calcular uma matriz de distâncias baseada na identidade, que 
    representa a proporção de diferenças entre as sequências. A partir desta matriz, a árvore filogenética é inferida 
    recorrendo ao método UPGMA e posteriormente visualizada de forma gráfica, permitindo a interpretação das relações 
    evolutivas entre as sequências analisadas.
    """
    
    alignment = AlignIO.read(aligned_fasta, "fasta")

    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    for clade in tree.get_nonterminals():
        clade.name = None

    fig, ax = plt.subplots(figsize=(16, 9))
    plt.rcParams.update({"font.size": 15})

    Phylo.draw(
        tree,
        axes=ax,
        label_func=lambda clade: id_to_organism.get(clade.name)
        if clade.is_terminal()
        else None,
        do_show=False
    )

    x_left, x_right = ax.get_xlim()
    ax.set_xlim(x_left, x_right * 1.35)

    plt.show()

    return tree

if __name__ == "__main__":

    xml_file = "blast_result_8676102.xml"
    query_fasta = "gene_angela.fasta"
    email = "teu_email@email.com"

    mafft_path = r"C:\Users\luisp\Downloads\mafft-7.526-win64-signed\mafft-win\mafft.bat"
    # O caminho para o executável do MAFFT deve ser adaptado ao sistema e à instalação local do utilizador que executar o código.

    fasta_hits, id_to_organism = msa_from_blastn_xml(
        xml_file=xml_file,
        query_fasta=query_fasta,
        email=email,
        n_hits=10
    )

    aligned = run_mafft(fasta_hits, mafft_path)

    tree = phylo_tree_nucl(aligned, id_to_organism)


