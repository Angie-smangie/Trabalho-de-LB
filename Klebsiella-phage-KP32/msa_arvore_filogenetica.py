from Bio.Blast import NCBIXML
from Bio import Entrez, SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def msa(xml_file, output_fasta="hits_proteina.fasta", n_hits=10, email="teu_email@exemplo.com"):

    """
        Esta função recebe como entrada um ficheiro XML resultante de uma pesquisa BLAST
    e extrai os N melhores alinhamentos (hits), obtendo os respetivos accession numbers.
    De seguida, utiliza esses accessions para comunicar com o NCBI através do Entrez,
    fazendo o download das sequências proteicas correspondentes em formato FASTA.
    O resultado final é um ficheiro FASTA contendo apenas as sequências homólogas
    mais relevantes, pronto a ser utilizado em análises posteriores como alinhamento
    múltiplo de sequências (MSA) e construção de árvores filogenéticas.
    """

    Entrez.email = email

    with open(xml_file) as f:
        blast_record = NCBIXML.read(f)

    accessions = [
        aln.accession
        for aln in blast_record.alignments[:n_hits]]

    print("Accessions extraídos:")
    for acc in accessions:
        print(" -", acc)

    handle = Entrez.efetch(
        db="protein",
        id=",".join(accessions),
        rettype="fasta",
        retmode="text")

    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()

    SeqIO.write(records, output_fasta, "fasta")

    print(f"\nFASTA criado: {output_fasta}")
    print(f"Total de sequências: {len(records)}")

    return output_fasta


def phylo_tree(aligned_fasta, matrix="blosum62", draw_ascii=True):

    """
        Esta função recebe como entrada um ficheiro FASTA contendo sequências já alinhadas
    por um método de alinhamento múltiplo, como o MAFFT. A partir desse alinhamento,
    é   calculada uma matriz de distâncias utilizando uma matriz de substituição adequada
    para proteínas (BLOSUM62). Com base nessa matriz de distâncias, é construída uma    
    árvore filogenética recorrendo ao método UPGMA, permitindo visualizar as relações
    evolutivas entre as sequências analisadas de forma hierárquica.
    """


    alignment = AlignIO.read(aligned_fasta, "fasta")

    calculator = DistanceCalculator(matrix)
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    if draw_ascii:
        Phylo.draw_ascii(tree)

    return tree