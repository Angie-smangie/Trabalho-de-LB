import time
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def buscar_genoma(email, cod_NCBI, api_key=None):
    """
    Descarrega um genoma do NCBI e guarda-o em formato GenBank (.gb)
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=cod_NCBI,
            rettype="gb",
            retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()

        with open(f"{cod_NCBI}.gb", "w") as output:
            SeqIO.write(record, output, "genbank")

        print(f"Genoma {cod_NCBI} descarregado com sucesso.")

        time.sleep(1)

    except Exception as e:
        print("Erro ao contactar o NCBI.")
        print(e)
    

def extrair_gene(cod_NCBI, gene_id):
    """
    Extrai informação completa de um gene (feature 'gene')
    com base no GeneID.

    É necessário utilizar o buscar_genoma antes ou obter o 
    ficheiro .gb antes de proceder com este.
    """
    with open(f"{cod_NCBI}.gb", "r") as ficheiro:
        genome = SeqIO.read(ficheiro, "genbank")

    for feature in genome.features:
        if feature.type == "gene":
            if gene_id in feature.qualifiers.get("db_xref", []):

                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand

                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                gene_name = feature.qualifiers.get("gene", [""])[0]
                db_xref = feature.qualifiers.get("db_xref", [])

                # sequência de DNA
                gene_seq = genome.seq[start:end]

                record = SeqRecord(
                    gene_seq,
                    id=gene_id,
                    description=f"locus_tag={locus_tag}"
                )

                with open(f"{gene_id}_gene.fasta", "w") as output:
                    SeqIO.write(record, output, "fasta")

                print("\n=== GENE ===")
                print(f"GeneID: {gene_id}")
                print(f"Locus tag: {locus_tag}")
                print(f"Gene name: {gene_name}")
                print(f"Localização: {start}..{end} (strand {strand})")
                print(f"db_xref: {db_xref}")

                return {
                    "GeneID": gene_id,
                    "locus_tag": locus_tag,
                    "gene_name": gene_name,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "db_xref": db_xref
                }

    print(f"Gene {gene_id} não encontrado.")
    return None

def extrair_cds(cod_NCBI, gene_id):
    """
    Extrai informação completa da CDS associada a um GeneID.
    
    É necessário utilizar o buscar_genoma antes ou obter o 
    ficheiro .gb antes de proceder com este.
    """
    with open(f"{cod_NCBI}.gb", "r") as ficheiro:
        genome = SeqIO.read(ficheiro, "genbank")

    for feature in genome.features:
        if feature.type == "CDS":
            if gene_id in feature.qualifiers.get("db_xref", []):

                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                product = feature.qualifiers.get("product", [""])[0]
                protein_id = feature.qualifiers.get("protein_id", [""])[0]
                translation = feature.qualifiers.get("translation", [""])[0]
                note = feature.qualifiers.get("note", [""])[0]
                ec = feature.qualifiers.get("EC_number", [""])[0]
                transl_table = feature.qualifiers.get("transl_table", [""])[0]

                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = feature.location.strand
                protein_len = len(translation)

                record = SeqRecord(
                    Seq(translation),
                    id=protein_id if protein_id else gene_id,
                    description=product
                )

                with open(f"{gene_id}_protein.fasta", "w") as output:
                    SeqIO.write(record, output, "fasta")

                print("\n=== CDS ===")
                print(f"GeneID: {gene_id}")
                print(f"Locus tag: {locus_tag}")
                print(f"Protein ID: {protein_id}")
                print(f"Produto: {product}")
                print(f"Comprimento: {protein_len} aa")
                print(f"Localização: {start}..{end} (strand {strand})")
                if ec:
                    print(f"EC number: {ec}")
                if note:
                    print(f"Nota: {note}")

                return {
                    "GeneID": gene_id,
                    "locus_tag": locus_tag,
                    "protein_id": protein_id,
                    "product": product,
                    "length_aa": protein_len,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "EC_number": ec,
                    "note": note,
                    "transl_table": transl_table
                }

    print(f"CDS associada a {gene_id} não encontrada.")
    return None
