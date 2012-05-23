from geneutils import *

DB_DIR = './'
MAIN_SPECIES='S288C_cerevisiae'
SUBJECT_DB = 'schizosaccharomyces_pombe_972h-_2'

# RUN                                                                                                                                                                                                       

# GENERATE BLAST DATABASES, GIVEN FASTA FILES                                                                                                                                                               
generate_blast_database(DB_DIR)

# SAVE QUERY AND SUBJECT DBS INTO A TWO-KEY DICTIONARY, SEARCHABLE BY SPECIES AND ORF NAME                                                                                                                  
generate_python_pep_database(DB_DIR, [SUBJECT_DB, MAIN_SPECIES])

def series():
    x=-10
    while x < -1:
        yield 10**x
        x += 1

for ev in series():
    counts=0
    blast('P', DB_DIR + MAIN_SPECIES + "_pep.fsa", DB_DIR + SUBJECT_DB, outname=MAIN_SPECIES + "-" + SUBJECT_DB + "---" + str(ev) + ".blastp.csv", evalue=ev)
    for (query_orf, subject_orf) in qseq_sseq_pairs(MAIN_SPECIES + '-' + SUBJECT_DB + "---" + str(ev) + '.blastp.csv'):
        sseq_fsa = LOCAL_PEP_DATABASE[SUBJECT_DB, subject_orf].format('fasta')
        if reverse_blast_check('P', MAIN_SPECIES, query_orf, sseq_fsa, evalue=ev):
            append_to_file(MAIN_SPECIES + '-' + SUBJECT_DB + "---" + str(ev) + '.blastp.csv.corrected', query_orf + "," + subject_orf + "\n")
            counts += 1
            print(counts)
    print("Finished with EVALUE = " + str(ev))
    append_to_file("STATS", str(ev) + "\t" + str(counts) + "\n")
