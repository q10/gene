import glob, shelve
from geneutils import *

DB_DIR = '../Proteomes/Candida_proteomes/'
BLASTP_RESULTS_DIR = 'resultsp/'

#MAIN_SPECIES = 'S288C'
#SUBJECT_DBS = ['AWRI1631_ABSV01000000', 'AWRI796_ADVS01000000', 'CBS7960_AEWL01000000', \
#'CLIB215_AEWP01000000', 'CLIB324_AEWM01000000', 'CLIB382_AFDG01000000', 'EC1118_PRJEA37863', \
#'EC9-8_AGSJ01000000', 'FL100_AEWO01000000', 'FostersB_AEHH01000000', 'FostersO_AEEZ01000000', \
#'Kyokai7_BABQ01000000', 'LalvinQA23_ADVV01000000', 'M22_ABPC01000000', 'PW5_AFDC01000000', \
#'RM11-1a_AAEG01000000', 'Sigma1278b_ACVY01000000', 'JAY291_ACFL01000000', 'T7_AFDE01000000', \
#'T73_AFDF01000000', 'UC5_AFDD01000000', 'Vin13_ADXC01000000', 'VL3_AEJS01000000', 'Y10_AEWK01000000', \
#'YJM269_AEWN01000000', 'YJM789_AAFW02000000', 'YPS163_ABPD01000000', 'W303_MPG_2011']

MAIN_SPECIES = 'S288C_cerevisiae__shift1'#, 'S288C_cerevisiae__shift2'
SUBJECT_DBS = ['candida_albicans_wo-1_1', 'candida_guilliermondii_1', 'candida_lusitaniae_1', \
'candida_parapsilosis_1', 'candida_tropicalis_3', 'debaryomyces_hansenii_1', 'lodderomyces_elongisporus_1']

#SUBJECT_DBS = ['K_waltii', 'L_kluyveri', 'S_bayanus', 'S_castellii', 'S_cerevisiae', 'S_kudriavzevii', \
#'S_mikatae', 'S_paradoxus', 'Scastellii_040406']

#SUBJECT_DBS = ['schizosaccharomyces_cryophilus_mito_3', 'schizosaccharomyces_cryophilus_oy26_3', \
#'schizosaccharomyces_japonicus_yfs275_4', 'schizosaccharomyces_japonicus_yfs275_mitochondria_1', \
#'schizosaccharomyces_octosporus_5', 'schizosaccharomyces_octosporus_mito', 'schizosaccharomyces_pombe_972h-_2', \
#'schizosaccharomyces_pombe_972h-_mitochondria']

# RUN

# GENERATE BLAST DATABASES, GIVEN FASTA FILES
generate_blast_database(DB_DIR)

# SAVE QUERY AND SUBJECT DBS INTO A TWO-KEY DICTIONARY, SEARCHABLE BY SPECIES AND ORF NAME
generate_python_pep_database(DB_DIR, SUBJECT_DBS + [MAIN_SPECIES])

'''
execute('mkdir -p ' + BLASTP_RESULTS_DIR + ' && rm -rf ' + BLASTP_RESULTS_DIR + '/*')
for entry in fasta_entries(DB_DIR + MAIN_SPECIES + '_pep.fsa'):
    write_file(BLASTP_RESULTS_DIR + entry.id, entry.format('fasta'))
    print("created pre-MUSCLE fasta file for " + entry.id)
'''
# RUN FULL BLASTP BETWEEN S288C AND EACH OF THE OTHER STRAINS
for DB_NAME in SUBJECT_DBS:
    blast('P', DB_DIR + MAIN_SPECIES + "_pep.fsa", DB_DIR + DB_NAME, outname=MAIN_SPECIES + "--" + DB_NAME + ".blastp.csv")

# TAKE BLASTP RESULTS, RUN REVERSE BLASTP, AND APPEND THE SEQUENCES TO THE APPROPRIATE FASTA FILES FOR MUSLCE.
# TAKES ONLY THE FIRST MATCH INTO THE PRE-MUSCLE FILE
for DB_NAME in SUBJECT_DBS:
    for (query_orf, subject_orfs_set) in qseq_sseq_sets(MAIN_SPECIES + '--' + DB_NAME + '.blastp.csv'):
        for subject_orf in subject_orfs_set:
            sseq_fsa = LOCAL_PEP_DATABASE[DB_NAME, subject_orf].format('fasta')
            if reverse_blast_check('P', DB_DIR + MAIN_SPECIES, query_orf, sseq_fsa):
                append_to_file(l2str(MAIN_SPECIES, '-', DB_NAME, "---", '.blastp__corrected.csv'), l2str(query_orf, ",", subject_orf, "\n"))
                break
    print("Finished processing " + DB_NAME + " matches")


titlerow = []
atable = []
for entry in fasta_entries(DB_DIR + MAIN_SPECIES + '_pep.fsa'):
    atable[entry.id] = []

count = 0
for DB_NAME in SUBJECT_DBS:
    count += 1
    titlerow.append(DB_NAME)
    for k, v in qseq_sseq_pairs(l2str(MAIN_SPECIES, '-', DB_NAME, "---", '.blastp__corrected.csv')):
        atable[k].append(v)
    for a in atable: # dangerous - could be one-to-many mapping
        if len(a) < count:
            a.append['']
        elif len(a) > count:
            a.pop()

append_to_file('COMBINED.CSV', ','.join(titlerow) + "\n")
for k, vtable in atable:
    append_to_file('COMBINED.CSV', l2str(k, ',', ','.join(vtable), "\n"))

# RUN MUSCLE
'''
failed_files = []
for multi_fasta_file in glob.glob(os.path.join(BLASTP_RESULTS_DIR + "*")):
    try:
        muscle(multi_fasta_file)
    except:
        failed_files.append(multi_fasta_file)

print('ALL DONE')
if failed_files:
    print('\n\nTHE FOLLOWING PRE-MUSCLE FILES FAILED:')
    for i in failed_files:
        print(i)
'''