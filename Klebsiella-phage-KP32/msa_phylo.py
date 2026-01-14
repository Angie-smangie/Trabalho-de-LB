from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import subprocess
import matplotlib.pyplot as plt


def msa_from_blastn_xml(
    xml_file: str,
    output_fasta: str = "hits_nucl.fasta",
    n_hits: int = 10,
    email: str = "teu_email@exemplo.com",
):
    """
    Os resultados BLASTn em formato XML foram processados de forma automática, sendo 
    selecionados os n melhores alinhamentos. Para cada hit, foi considerado apenas o 
    HSP com maior score, correspondente à região de maior similaridade entre a sequência
    query e a sequência da base de dados. Com base nas coordenadas desse alinhamento, foi
    extraída exclusivamente a região homóloga de cada sequência alvo, recorrendo à base 
    de dados nuccore do NCBI. As regiões obtidas foram reunidas num ficheiro FASTA multi-sequência,
    que serviu de base para o alinhamento múltiplo e para a análise filogenética subsequente
    """

    Entrez.email = email

    with open(xml_file) as f:
        blast_record = NCBIXML.read(f)

    if not blast_record.alignments:
        raise ValueError("O XML não tem alinhamentos BLAST (alignments vazios).")

    records = []
    print("Hits extraídos (com coordenadas do 1º HSP):")

    for aln in blast_record.alignments[:n_hits]:
        acc = aln.accession

        if not aln.hsps:
            print(f" - {acc}: sem HSPs (ignorado)")
            continue

        hsp = aln.hsps[0]
        s_start = int(hsp.sbjct_start)
        s_end = int(hsp.sbjct_end)

        start = min(s_start, s_end)
        end = max(s_start, s_end)

        strand = 1 if s_start <= s_end else 2

        print(f" - {acc}: {start}-{end} (strand={'+' if strand == 1 else '-'})")

        handle = Entrez.efetch(
            db="nuccore",
            id=acc,
            rettype="fasta",
            retmode="text",
            seq_start=start,
            seq_stop=end,
            strand=strand,
        )

        rec = SeqIO.read(handle, "fasta")
        handle.close()

        rec.id = acc
        rec.name = acc
        rec.description = f"{acc} region {start}-{end} strand {strand}"

        records.append(rec)

    if not records:
        raise ValueError("Não foi possível obter nenhuma sequência do NCBI (records vazio).")

    SeqIO.write(records, output_fasta, "fasta")
    print(f"\nFASTA criado: {output_fasta}")
    print(f"Total de sequências: {len(records)}")

    return output_fasta


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


def phylo_tree_nucl(aligned_fasta: str):
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

    Phylo.draw(tree)
    plt.show()

    return tree


if __name__ == "__main__":

    xml_file = ".xml" 
    email = "teu_email@email.com"

    mafft_path = r"C:\Users\utilizador..."
    # O caminho para o executável do MAFFT deve ser adaptado ao sistema e à instalação local do utilizador que executar o código.

    fasta_hits = msa_from_blastn_xml(xml_file, email=email, n_hits=10)
    aligned = run_mafft(fasta_hits, mafft_path)
    tree = phylo_tree_nucl(aligned)from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import subprocess
import matplotlib.pyplot as plt


def msa_from_blastn_xml(
    xml_file: str,
    output_fasta: str = "hits_nucl.fasta",
    n_hits: int = 10,
    email: str = "teu_email@exemplo.com",
):
    """
    Os resultados BLASTn em formato XML foram processados de forma automática, sendo 
    selecionados os n melhores alinhamentos. Para cada hit, foi considerado apenas o 
    HSP com maior score, correspondente à região de maior similaridade entre a sequência
    query e a sequência da base de dados. Com base nas coordenadas desse alinhamento, foi
    extraída exclusivamente a região homóloga de cada sequência alvo, recorrendo à base 
    de dados nuccore do NCBI. As regiões obtidas foram reunidas num ficheiro FASTA multi-sequência,
    que serviu de base para o alinhamento múltiplo e para a análise filogenética subsequente
    """

    Entrez.email = email

    with open(xml_file) as f:
        blast_record = NCBIXML.read(f)

    if not blast_record.alignments:
        raise ValueError("O XML não tem alinhamentos BLAST (alignments vazios).")

    records = []
    print("Hits extraídos (com coordenadas do 1º HSP):")

    for aln in blast_record.alignments[:n_hits]:
        acc = aln.accession

        if not aln.hsps:
            print(f" - {acc}: sem HSPs (ignorado)")
            continue

        hsp = aln.hsps[0]
        s_start = int(hsp.sbjct_start)
        s_end = int(hsp.sbjct_end)

        start = min(s_start, s_end)
        end = max(s_start, s_end)

        strand = 1 if s_start <= s_end else 2

        print(f" - {acc}: {start}-{end} (strand={'+' if strand == 1 else '-'})")

        handle = Entrez.efetch(
            db="nuccore",
            id=acc,
            rettype="fasta",
            retmode="text",
            seq_start=start,
            seq_stop=end,
            strand=strand,
        )

        rec = SeqIO.read(handle, "fasta")
        handle.close()

        rec.id = acc
        rec.name = acc
        rec.description = f"{acc} region {start}-{end} strand {strand}"

        records.append(rec)

    if not records:
        raise ValueError("Não foi possível obter nenhuma sequência do NCBI (records vazio).")

    SeqIO.write(records, output_fasta, "fasta")
    print(f"\nFASTA criado: {output_fasta}")
    print(f"Total de sequências: {len(records)}")

    return output_fasta


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


def phylo_tree_nucl(aligned_fasta: str):
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

    Phylo.draw(tree)
    plt.show()

    return tree


if __name__ == "__main__":

    xml_file = "blastn.xml" 
    email = "teu_email@email.com"

    mafft_path = r"C:\Users\luisp\Downloads\mafft-7.526-win64-signed\mafft-win\mafft.bat"
    # O caminho para o executável do MAFFT deve ser adaptado ao sistema e à instalação local do utilizador que executar o código.

    fasta_hits = msa_from_blastn_xml(xml_file, email=email, n_hits=10)
    aligned = run_mafft(fasta_hits, mafft_path)
    tree = phylo_tree_nucl(aligned)




