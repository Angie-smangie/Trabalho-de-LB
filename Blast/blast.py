#é necessário que biopython esteja instalado no pc
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from io import StringIO

def blast_aa(*seqs, program="blastp", database="nr", hitlist_size=10, save_xml=True):
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

