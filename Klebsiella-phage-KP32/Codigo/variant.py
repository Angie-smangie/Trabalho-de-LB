def buscar_gene(cod_NCBI,gene):
    """
    Esta função permite ir buscar um gene a um genoma .gb, e devolve informações sobre o mesmo, juntamente com um ficheiro fasta com a sequência do mesmo.

    cod_NCBI: aceita um código de genoma do NCBI
    gene: aceita um GeneID no formato 'GeneID:XXXXXXX'
    """
    with open(f"{cod_NCBI}.gb", "r") as ficheiro:
        genoma = SeqIO.read(ficheiro, "genbank")
    for i, feature in enumerate(genoma.features):
        if feature.type == "gene":
            db_xrefs = feature.qualifiers.get("db_xref", [])
            if gene in db_xrefs:
                cod=re.search(r"[\d]+",gene)
                start = int(feature.location.start)
                end = int(feature.location.end)
                location = genoma.seq[start:end+1]
                locus = feature.qualifiers.get("locus_tag", [""])[0]
                record = SeqRecord(seq=location,id=gene,description=locus)
                with open(f"{cod[0]}.fasta","w") as output:
                    SeqIO.write(record, output, "fasta")
                print(f"{gene} encontrado")
                print("Posição no ficheiro:", i)
                print("Localização:", feature.location)
                return feature

    print(f"{gene} não encontrado.")
    return None

def buscar_traducao(cod_NCBI,gene):
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
                cod=re.search(r"[\d]+",gene)
                locus = feature.qualifiers.get("locus_tag", [""])[0]
                protein_seq = feature.qualifiers.get("translation", [""])[0]
                record = SeqRecord(seq=Seq(protein_seq),id=gene,description=locus)
                with open(f"{cod[0]}_translated.fasta","w") as output:
                    SeqIO.write(record, output, "fasta")
                print(f"CDS de {gene} encontrado")
                print("Posição no ficheiro:", i)
                print("Localização:", feature.location)
                return feature
            
def mostrar_resultados(registo_blast):
    """
    Processa e exibe no ecrã um resumo estatístico dos melhores hits de um registo BLAST.
    
    A função calcula a percentagem de identidade e de cobertura da query para cada hit,
    imprimindo o Accession, descrição, E-value e as métricas calculadas.

    Parâmetros:
    - registo_blast: objeto Blast do Biopython (Bio.Blast.Record.Blast) contendo
      os resultados de uma query individual.

    Retorna:
    - None: Apenas imprime os dados formatados na consola.
    """
    
    # Extrair o objeto se for passado um dicionário
    if isinstance(registo_blast, dict):
        if not registo_blast:
            print("Objeto vazio.")
            return
        # Pega no primeiro valor do dicionário
        registo_blast = list(registo_blast.values())[0]

    # Verifica se temos resultados válidos
    if not registo_blast:
        print("Objeto vazio ou inválido.")
        return

    meus_hits = get_top_hits(registo_blast, top_n=5)
    
    # Tamanho da Query
    tamanho_total_query = registo_blast.query_letters

    # --- 2. CORREÇÃO AQUI ---
    # Usamos apenas .query ou uma string fixa caso falhe. Removemos .id
    nome_query = getattr(registo_blast, 'query', 'Query Desconhecida')
    # ------------------------

    print(f"--- ANÁLISE PARA: {nome_query} ---")

    for hit in meus_hits:
        hsp = hit['hsps'][0]
        
        len_alinhamento = hsp.get('aln_span', hsp.get('align_length'))
        num_identidades = hsp['identities']
        
        identidade_pc = (num_identidades / len_alinhamento) * 100
    
        if tamanho_total_query > 0:
            cobertura_pc = (len_alinhamento / tamanho_total_query) * 100
        else:
            cobertura_pc = 0
            
        if cobertura_pc > 100: cobertura_pc = 100.0

        print(f"HIT: {hit['accession']}")
        descricao = hit.get('description', hit.get('definition', 'Sem descrição'))
        print(f"Descrição: {descricao[:80]}...") 
        print(f"  > E-value:    {hsp['e_value']:.2e}")   
        print(f"  > Identidade: {identidade_pc:.2f}%")   
        print(f"  > Cobertura:  {cobertura_pc:.2f}%")
        print("-" * 40)

def guardar_para_mega(resultado_blast, email, database, nome_do_ficheiro, quantidade, ficheiro_original):
    """
    Processa os resultados do BLAST para extrair variantes genéticas e guardar num ficheiro FASTA, 
    preparando os dados para alinhamento no MEGA.

    A função realiza os seguintes passos:
    1. Se fornecido, adiciona a sequência original (query) no início do ficheiro para servir de referência.
    2. Percorre os hits do BLAST e ignora sequências 100% idênticas à query.
    3. Descarrega as sequências das variantes encontradas:
        - Se database='protein': baixa as sequências de aminoácidos completas.
        - Se database='nucleotide': baixa apenas a região alinhada (o gene), respeitando as coordenadas.

    Parâmetros:
    - resultado_blast: Objeto com o resultado de uma query BLAST.
    - email: O email (obrigatório para o NCBI).
    - database: "protein" ou "nucleotide".
    - nome_do_ficheiro: Nome do ficheiro de saída (.fasta).
    - quantidade: Número máximo de variantes a guardar.
    - ficheiro_original: Caminho para o ficheiro FASTA com a sequência query original.
    """
    Entrez.email = email
    
    lista_ids_proteinas = []
    
    with open(nome_do_ficheiro, "w") as ficheiro_saida:
        
        # Ficheiro original
        try:
            with open(ficheiro_original, "r") as f_orig:
                conteudo = f_orig.read().strip()
                
                if conteudo:
                    ficheiro_saida.write(conteudo)
                    ficheiro_saida.write("\n\n") 
            
            print("A query original foi adicionada com sucesso.")
        except Exception as e:
            print(f"Não consegui ler o ficheiro original: {e}")
        
        conta = 0
        
        if isinstance(resultado_blast, dict):
             if not resultado_blast: return
             resultado_blast = list(resultado_blast.values())[0]

        for alinhamento in resultado_blast.alignments:
            if conta >= quantidade:
                break
            
            hsp = alinhamento.hsps[0]
            
            # Ignorar 100% iguais (Variant Calling)
            if hsp.identities == hsp.align_length:
                print(f"-> A ignorar {alinhamento.accession} (é 100% igual)")
                continue 

            accession = alinhamento.accession
            
            if database == "nucleotide":
                # Para nucleótidos, baixamos região específica
                start = hsp.sbjct_start
                end = hsp.sbjct_end
                
                if start < end:
                    s_start, s_end, strand_val = start, end, 1
                else:
                    s_start, s_end, strand_val = end, start, 2
                
                try:
                    handle = Entrez.efetch(
                        db="nucleotide", 
                        id=accession, 
                        rettype="fasta", 
                        retmode="text",
                        seq_start=s_start, 
                        seq_stop=s_end,
                        strand=strand_val
                    )
                    
                    # Lê e limpa espaços extra
                    seq_dados = handle.read().strip()
                    
                    # Escreve a sequência e força linha em branco
                    ficheiro_saida.write(seq_dados)
                    ficheiro_saida.write("\n\n")
                    
                    handle.close()
                    conta += 1
                except Exception as e:
                    print(f"Erro ao baixar região de {accession}: {e}")

            else:
                lista_ids_proteinas.append(accession)
                conta += 1

        # Processamento de Proteínas
        if database == "protein" and lista_ids_proteinas:
            try:
                handle = Entrez.efetch(
                    db="protein", 
                    id=lista_ids_proteinas, 
                    rettype="fasta", 
                    retmode="text"
                )
                
                dados_prot = handle.read().strip()
                ficheiro_saida.write(dados_prot)
                ficheiro_saida.write("\n\n")
                
                handle.close()
                print("Sucesso! Proteínas guardadas.")
            except Exception as e:
                print(f"Erro ao baixar proteínas: {e}")
        
        elif database == "nucleotide" and conta > 0:
            print("Regiões de genes guardadas.")
        
        elif conta == 0:
            print("Aviso: Nenhuma sequência foi guardada.")
