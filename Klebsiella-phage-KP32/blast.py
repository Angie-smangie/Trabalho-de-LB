#é necessário que biopython esteja instalado no pc
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from io import StringIO

def blast_aa(*seqs, program="blastp", database="nr", hitlist_size=10, save_xml=True):

    """
    Esta função realiza um BLAST de sequências de aminoácidos no NCBI
    e guarda os resultados num ficheiro XML, além de retornar os registros
    processados pelo Biopython.

    seqs : aceita uma ou mais sequências de aminoácidos, que podem ser strings
           ou objetos SeqRecord do Biopython.
    program : string opcional, programa BLAST a utilizar (padrão: "blastn").
    database : string opcional, base de dados do NCBI a pesquisar (padrão: "nt").
    hitlist_size : inteiro opcional, número máximo de hits a recuperar por sequência (padrão: 10).
    save_xml : booleano opcional, se True guarda o resultado BLAST em ficheiro XML (padrão: True).

    A função imprime no ecrã os 10 melhores alinhamentos da última sequência
    processada, incluindo Accession, Definition e E-value, e devolve um dicionário
    com os SeqIDs como chave e objetos BlastRecord como valor.
    """
    
    resultados = {}
    for i, s in enumerate(seqs, start=1):
        if isinstance(s, SeqRecord):
            seq_str = str(s.seq)
            seq_id = s.id
        elif isinstance(s, str):
            seq_str = s
            seq_id = f"seq_{i}"
        else:
            raise TypeError("A sequência tem de ser string ou SeqRecord")

        print(f"BLAST para {seq_id}...")
        handle = NCBIWWW.qblast(program, database, seq_str, hitlist_size=hitlist_size)

        xml_text = handle.read()

        if save_xml:
            xml_file = f"{seq_id}.xml"
            with open(xml_file, "w", encoding="utf-8") as out:
                out.write(xml_text)
            print(f"XML criado: {xml_file}")

        blast_record = NCBIXML.read(StringIO(xml_text))
        resultados[seq_id] = blast_record
    
    for i in range(10):
        alignment = blast_record.alignments[i]
        print ("Accession: " + alignment.accession)
        print ("Definition: " + alignment.hit_def)
        for hsp in alignment.hsps:
            print ("E-value: ", hsp.expect)
        
    return resultados
   

def blast_nucl(*seqs, program="blastn", database="nt", hitlist_size=10, save_xml=True):
   
    """
    Esta função realiza um BLAST de sequências de nucleótidos no NCBI
    e guarda os resultados num ficheiro XML, além de retornar os registros
    processados pelo Biopython.

    seqs : aceita uma ou mais sequências de nucleótidos, que podem ser strings
           ou objetos SeqRecord do Biopython.
    program : string opcional, programa BLAST a utilizar (padrão: "blastn").
    database : string opcional, base de dados do NCBI a pesquisar (padrão: "nt").
    hitlist_size : inteiro opcional, número máximo de hits a recuperar por sequência (padrão: 10).
    save_xml : booleano opcional, se True guarda o resultado BLAST em ficheiro XML (padrão: True).

    A função imprime no ecrã os 10 melhores alinhamentos da última sequência
    processada, incluindo Accession, Definition e E-value, e devolve um dicionário
    com os SeqIDs como chave e objetos BlastRecord como valor.
    """
    
    resultados = {}
    for i, s in enumerate(seqs, start=1):
        if isinstance(s, SeqRecord):
            seq_str = str(s.seq)
            seq_id = s.id
        elif isinstance(s, str):
            seq_str = s
            seq_id = f"seq_{i}"
        else:
            raise TypeError("A sequência tem de ser string ou SeqRecord")

        print(f"BLAST para {seq_id}...")
        handle = NCBIWWW.qblast(program, database, seq_str, hitlist_size=hitlist_size)

        xml_text = handle.read()

        if save_xml:
            xml_file = f"{seq_id}.xml"
            with open(xml_file, "w", encoding="utf-8") as out:
                out.write(xml_text)
            print(f"XML criado: {xml_file}")

        blast_record = NCBIXML.read(StringIO(xml_text))
        resultados[seq_id] = blast_record

    for i in range(10):
        alignment = blast_record.alignments[i]
        print ("Accession: " + alignment.accession)
        print ("Definition: " + alignment.hit_def)
        for hsp in alignment.hsps:
            print ("E-value: ", hsp.expect)

    return resultados


