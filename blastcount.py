from geneutils import *

DB_DIR = '../Species_Genomes/Ascomycota/'
MAIN_SPECIES='S288C_cerevisiae'
SUBJECT_DB = 'schizosaccharomyces_pombe_972h-_2'

# RUN

# GENERATE BLAST DATABASES, GIVEN FASTA FILES
generate_blast_database()

# SAVE QUERY AND SUBJECT DBS INTO A TWO-KEY DICTIONARY, SEARCHABLE BY SPECIES AND ORF NAME
generate_python_cds_database([SUBJECT_DB, MAIN_SPECIES])

def series():
    x=-10
    while x < -1:
        yield 10**x
        x += 1

for ev in series():
    counts=0
    blastn(DB_DIR + MAIN_SPECIES + "_cds.fsa", DB_DIR + SUBJECT_DB, outname=MAIN_SPECIES + "-" + SUBJECT_DB + "---" + str(ev) + ".blastn.csv", evalue=ev)
    for (query_orf, subject_orf) in qseq_sseq_pairs(MAIN_SPECIES + '-' + SUBJECT_DB + "---" + str(ev) + '.blastn.csv'):
        sseq_fsa = LOCAL_CDS_DATABASE[SUBJECT_DB, subject_orf].format('fasta')
        if reverse_blastn_check(MAIN_SPECIES, query_orf, sseq_fsa):
            append_to_file(MAIN_SPECIES + '-' + SUBJECT_DB + "---" + str(ev) + '.blastn.csv.corrected', sseq_fsa)
            counts += 1
    print("Finished with EVALUE = " + str(ev))
    append_to_file("STATS", str(ev) + "\t" + str(counts) + "\n")
