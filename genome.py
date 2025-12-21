#É necessário alterar o código de maneira a ir buscar as anotações e toda a info relevante

import Bio
from Bio import Entrez
from Bio import SeqIO

def buscar_genoma(email,cod_NCBI):
    """
    Esta função permite ir buscar um genoma ao NCBI, que guarda num ficheiro .gb.
    
    email: aceita um email em formato string
    cod_NCBI: aceita um código de genoma do NCBI em formato string
    """
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=cod_NCBI, rettype="gb", retmode="text")
    result = SeqIO.read(handle, "genbank")
    handle.close()
    with open(f"{cod_NCBI}.gb","w") as genome:
        SeqIO.write(result,genome,"genbank")
    

def buscar_gene(cod_NCBI,gene):
    """
    Esta função permite ir buscar um gene a um genoma .gb, e devolve informações sobre o mesmo.

    cod_NCBI: aceita um código de genoma do NCBI
    gene: aceita um GeneID no formato 'GeneID:XXXXXXX'
    """
    with open(f"{cod_NCBI}.gb", "r") as ficheiro:
        genome = SeqIO.read(ficheiro, "genbank")

    for i, feature in enumerate(genome.features):
        if feature.type == "gene":
            db_xrefs = feature.qualifiers.get("db_xref", [])
            if gene in db_xrefs:
                print(f"{gene} encontrado")
                print("Posição no ficheiro:", i)
                print("Localização:", feature.location)
                return feature

    print(f"{gene} não encontrado.")
    return None

def buscar_CDS(cod_NCBI,gene):
    """
    Esta função permite ir buscar a CDS associada a um gene de um genoma .gb, e devolve informações sobre o mesmo.

    cod_NCBI: aceita um código de genoma do NCBI
    gene: aceita um GeneID no formato 'GeneID:XXXXXXX'
    """
    with open(f"{cod_NCBI}.gb", "r") as ficheiro:
        genome = SeqIO.read(ficheiro, "genbank")

    for i, feature in enumerate(genome.features):
        if feature.type == "CDS":
            db_xrefs = feature.qualifiers.get("db_xref", [])
            if gene in db_xrefs:
                print(f"CDS de {gene} encontrado")
                print("Posição no ficheiro:", i)
                print("Localização:", feature.location)
                return feature

    print(f"CDS para {gene} não encontrado.")
    return None