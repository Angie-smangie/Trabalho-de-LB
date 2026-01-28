from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from io import StringIO
import os
import re

def run_blast(input_data, program, database, hitlist_size=500, e_value=10.0, low_complexity_filter=False, save_xml=True):
    """
    Executa BLAST via NCBIWWW.qblast com parâmetros ajustados para resultados próximos ao site NCBI.

    Parâmetros:
    - input_data: str ou SeqRecord ou arquivo FASTA
    - program: 'blastn', 'blastp', etc.
    - database: 'nt', 'nr', etc.
    - hitlist_size: número máximo de hits a retornar (default 500)
    - e_value: cutoff do E-value (default 10.0)
    - low_complexity_filter: aplicar filtro de regiões de baixa complexidade (default False)
    - save_xml: salvar XML localmente

    Retorna:
    - dict: Dicionário onde cada chave é "seq_X" (ex: "seq_1", "seq_2") e cada valor é
      um objeto Blast record do Biopython (NCBIXML.Blast) correspondente à sequência processada.
      Esses objetos podem ser usados para extrair alinhamentos, HSPs e outras informações do BLAST.
    """
    
    # confirma se o input é correto
    if isinstance(input_data, str) and os.path.isfile(input_data):
        seq_records = list(SeqIO.parse(input_data, "fasta"))
    elif isinstance(input_data, str):
        seq_records = [SeqRecord(Seq(input_data), id="input_sequence")]
    elif isinstance(input_data, SeqRecord):
        seq_records = [input_data]
    else:
        raise TypeError("input_data inválido")

    resultados = {}

    for i, seq in enumerate(seq_records, start=1):
        seq_str = str(seq.seq)

        #nome do ficheiro para guardar
        match = re.search(r'\[GeneID=(\d+)\]', seq.description)
        if match:
                safe_id = match.group(1)
        else:
            safe_id = re.sub(r'[^A-Za-z0-9_\-]', '_', seq.id)
        xml_name = f"blast_result_{safe_id}_{program}.xml"



        print(f"A correr {program} -> {xml_name}")
        # faz o blast
        if program == "blastn":
            handle = NCBIWWW.qblast(
                program=program,
                database=database,
                sequence=seq_str,
                hitlist_size=hitlist_size,
                expect=e_value,
                entrez_query='NOT Klebsiella[All Fields] AND NOT phage[All Fields]',
                filter='L' if low_complexity_filter else 'F',  # 'F' desativa filtro
                format_type='XML'
            )
        else:
            handle = NCBIWWW.qblast(
                program=program,
                database=database,
                sequence=seq_str,
                hitlist_size=hitlist_size,
                expect=e_value,
                filter='L' if low_complexity_filter else 'F',  # 'F' desativa filtro
                format_type='XML'
            )
        xml_text = handle.read()
        handle.close()

         # salva o XML localmente(se na funçao tiver save_xml=True)
        if save_xml:
            with open(xml_name, "w", encoding="utf-8") as f:
                f.write(xml_text)
            
        #cria um objeto Blast record do biopython que é armazenado no dicionário de resultados
        blast_record = NCBIXML.read(StringIO(xml_text))
        resultados[f"seq_{i}"] = blast_record

    return resultados


def get_top_hits(blast_record, top_n=10):
    """
    Retorna os detalhes dos top N hits de um blast_record como uma lista de dicionários.

    Parâmetros:
    - blast_record: objeto BLAST do Biopython (NCBIXML.read)
    - top_n: número de hits a retornar (default 10)

    Retorna:
    - List[Dict]: cada elemento é um hit com os dados de seus HSPs
    """
    top_hits = []
    alignments = blast_record.alignments[:top_n]
    # ve os hits e aponta as chaves desse hit
    for i, aln in enumerate(alignments, start=1):
        hit_data = {
            "hit_number": i,
            "accession": aln.accession,
            "definition": aln.hit_def,
            "hsps": []
        }
        # ve o hsp e atualiza as chaves desse hsp
        for hsp in aln.hsps:
            hsp_data = {
                "e_value": hsp.expect,
                "score": hsp.score,
                "bit_score": hsp.bits,
                "query_start": hsp.query_start,
                "query_end": hsp.query_end,
                "subject_start": hsp.sbjct_start,
                "subject_end": hsp.sbjct_end,
                "align_length": hsp.align_length,
                "identities": hsp.identities,
                "gaps": hsp.gaps
            }
            hit_data["hsps"].append(hsp_data)

        top_hits.append(hit_data)

    return top_hits
