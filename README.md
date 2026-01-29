# **Laboratórios de Bioinformática – Grupo 1**

Este repositório contém o código, dados e resultados desenvolvidos no âmbito do trabalho do **Grupo 1** da unidade curricular **Laboratórios de Bioinformática**, constituído por **Ângela Sousa**, **António Araújo**, **Filipe Kuhlmann** e **Luís Pedrosa**.

O mesmo tem como objetivo a análise de quatro genes de *Klebsiella phage* KP32, nomeadamente **KP32gp15**, **KP32gp31**, **KP32gp31** & **KP32gp38**.  O *K. phage* KP32 é um bacteriófago que infeta *Klebsiella pneumoniae*, uma bactéria **Gram-negativa**, **multirresistente** e frequentemente **hipervirulenta**.

*K. pneumoniae* tipicamente existe inofensivamente no sistema digestivo humano, mais concretamente no trato gastrointestinal. É, no entanto, um patógeno oportunista capaz de causar infeções graves como infeções do trato urinário, pneumonia, abcessos hepáticos e septimicia, especialmente em ambientes hospitalares. A emergência de estirpes multirresistentes torna o estudo de bacteriófagos particularmente relevante como potencial alternativa terapêutica.

## **Estrutura do Repositório**

### Pasta **Klebsiella-phage-KP32**

Contém os *scripts* principais desenvolvidos em Python, bem como dados de referência:

* codigo/
    * genome.py – aquisição de genomas .gb e aquisição sobre sobre genes pertencentes ao mesmo
    * blast.py – execução de BLAST
    * msa_phylo.py – alinhamento múltiplo de sequências e construção de árvores filogenéticas
    * exemplo_utiliz_cod.ipynb – notebook demonstrativo com exemplos de utilização do código desenvolvido
      
* NC_013647.1.gb – ficheiro GenBank correspondente ao genoma do *K. phage* KP32
* **Resultados da ferramenta PhagePromoter**, utilizados para identificação de promotores no genoma do bacteriófago

### Pastas dos Genes (**Gene_KP32gp15**, **Gene_KP32gp31**, **Gene_KP32gp37**, **Gene_KP32gp38**)

* XXX.fasta – sequência nucleotídica do gene
* XXX_translated.fasta – sequência de aminoácidos
* blast_result_kp32gpXX.xml e blast_result_kp32gpXX_blastp.xml - resultados BLASTN e BLASTP, respetivamente. Correspondem a análises com exclusão de *hits* de *K. phage*.
* hits_nucl_aligned_kp32gp38.fasta - resultados BLAST alinhados com MAFFT para construir a árvore filogenética, demonstrada em:
* phylo_kp32gpXX.jpeg
* sequencias_mega_proteinaXX - resultados BLASTP sem exclusão de *hits* de *K. phage*. Este ficheiro foi utilizado para se realizar o alinhamento MUSCLE no MEGA, disponível em:
* kp32gpXX_muscle.fas
* cdd_XX.png & cdd_XX_X.png - resultados do Conserved Domain Search do NCBI.

## **Análise dos Resultados**

A análise e interpretação detalhada dos resultados obtidos neste trabalho encontra-se disponível no seguinte link: https://padlet.com/asgdes/grupo-1-trabalho-de-laborat-rios-de-bioinform-tica-9p335t6w2qe4pp6j
