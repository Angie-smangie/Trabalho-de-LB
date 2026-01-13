from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import subprocess


def msa(xml_file, output_fasta="hits_proteina.fasta", n_hits=10, email="teu_email@exemplo.com"):
    """
    A função msa tem como objetivo processar um ficheiro XML resultante de uma pesquisa BLAST 
    e obter as sequências proteicas homólogas mais relevantes. A partir do ficheiro XML, a função 
    identifica os melhores alinhamentos e extrai os respetivos accession numbers. De seguida, 
    estabelece comunicação com a base de dados do NCBI através do módulo Entrez, permitindo o 
    descarregamento automático das sequências proteicas correspondentes em formato FASTA. O resultado 
    desta função é um ficheiro FASTA que reúne apenas as sequências homólogas selecionadas, constituindo 
    a base para análises subsequentes, nomeadamente o alinhamento múltiplo de sequências e a construção 
    de árvores filogenéticas.
    """

    Entrez.email = email

    with open(xml_file) as f:
        blast_record = NCBIXML.read(f)

    accessions = [aln.accession for aln in blast_record.alignments[:n_hits]]

    print("Accessions extraídos:")
    for acc in accessions:
        print(" -", acc)

    handle = Entrez.efetch(
        db="protein",
        id=",".join(accessions),
        rettype="fasta",
        retmode="text"
    )

    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    SeqIO.write(records, output_fasta, "fasta")

    print(f"\nFASTA criado: {output_fasta}")
    print(f"Total de sequências: {len(records)}")

    return output_fasta


def run_mafft(input_fasta, mafft_path, output_fasta="hits_proteina_aligned.fasta"):
    """
    A função run_mafft é responsável pela realização do alinhamento múltiplo de sequências (MSA), utilizando 
    o software MAFFT. Esta função recebe como entrada o ficheiro FASTA previamente gerado pela função msa, bem 
    como o caminho para o executável do MAFFT no sistema. O alinhamento é executado recorrendo à opção --auto, 
    que permite ao MAFFT selecionar automaticamente a estratégia de alinhamento mais adequada às sequências fornecidas. 
    O output desta etapa é um novo ficheiro FASTA contendo as sequências alinhadas, o qual será posteriormente utilizado 
    para a análise filogenética.
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


def phylo_tree(aligned_fasta, matrix="blosum62"):
    """
    A função phylo_tree tem como finalidade a construção de uma árvore filogenética a partir de um alinhamento 
    múltiplo de sequências proteicas. O ficheiro FASTA alinhado é utilizado para calcular uma matriz de distâncias 
    entre as sequências, com base numa matriz de substituição apropriada para proteínas (BLOSUM62). A partir desta 
    matriz de distâncias, é construída uma árvore filogenética recorrendo ao método UPGMA, permitindo inferir as 
    relações evolutivas entre as sequências analisadas. A árvore é visualizada graficamente no ambiente, facilitando 
    a interpretação das relações filogenéticas entre as sequências estudadas
    """

    alignment = AlignIO.read(aligned_fasta, "fasta")

    calculator = DistanceCalculator(matrix)
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    Phylo.draw(tree)
    plt.show()
    return tree


if __name__ == "__main__":

    xml_file = "file.xml"
    email = "teu_email@email.com"

    mafft_path = r"C:\Users\utilizador\..."
    # O caminho para o executável do MAFFT deve ser adaptado ao sistema e à instalação local do utilizador.

    fasta = msa(xml_file, email=email)
    aligned = run_mafft(fasta, mafft_path)
    tree = phylo_tree(aligned)



