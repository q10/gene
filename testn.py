import glob, shelve
from geneutils import *

BLASTN_RESULTS_DIR = 'resultsn/'
MAIN_SPECIES = 'S288C'
SUBJECT_DBS = ['AWRI1631_ABSV01000000', 'AWRI796_ADVS01000000', 'CBS7960_AEWL01000000', \
'CLIB215_AEWP01000000', 'CLIB324_AEWM01000000', 'CLIB382_AFDG01000000', 'EC1118_PRJEA37863', \
'EC9-8_AGSJ01000000', 'FL100_AEWO01000000', 'FostersB_AEHH01000000', 'FostersO_AEEZ01000000', \
'Kyokai7_BABQ01000000', 'LalvinQA23_ADVV01000000', 'M22_ABPC01000000', 'PW5_AFDC01000000', \
'RM11-1a_AAEG01000000', 'Sigma1278b_ACVY01000000', 'JAY291_ACFL01000000', 'T7_AFDE01000000', \
'T73_AFDF01000000', 'UC5_AFDD01000000', 'Vin13_ADXC01000000', 'VL3_AEJS01000000', 'Y10_AEWK01000000', \
'YJM269_AEWN01000000', 'YJM789_AAFW02000000', 'YPS163_ABPD01000000', 'W303_MPG_2011']

# RUN
'''
# GENERATE BLAST DATABASES, GIVEN FASTA FILES
generate_blast_database()
'''

'''
# SAVE QUERY AND SUBJECT DBS INTO A TWO-KEY DICTIONARY, SEARCHABLE BY SPECIES AND ORF NAME
generate_python_databases(SUBJECT_DBS + [MAIN_SPECIES])


load_python_databases()

execute('mkdir -p ' + BLASTN_RESULTS_DIR + ' && rm -rf ' + BLASTN_RESULTS_DIR + '*')
for entry in fasta_entries(DB_DIR + MAIN_SPECIES + '_cds.fsa'):
	write_file(BLASTN_RESULTS_DIR + entry.id, entry.format('fasta'))
	print("created pre-MUSCLE fasta file for " + entry.id)

# RUN FULL BLASTN BETWEEN S288C AND EACH OF THE OTHER STRAINS
for DB_NAME in SUBJECT_DBS:
    blastn(DB_DIR + MAIN_SPECIES + "_cds.fsa", DB_DIR + DB_NAME, outname=MAIN_SPECIES + "-" + DB_NAME + ".blastn.csv")

# TAKE BLASTP RESULTS, RUN REVERSE BLASTP, AND APPEND THE SEQUENCES TO THE APPROPRIATE FASTA FILES FOR MUSLCE
for DB_NAME in SUBJECT_DBS:
    for (query_orf, subject_orfs_set) in qseq_sseq_sets(MAIN_SPECIES + '-' + DB_NAME + '.blastn.csv'):
        for subject_orf in subject_orfs_set:
            sseq_fsa = LOCAL_CDS_DATABASE[DB_NAME, subject_orf].format('fasta')
            if reverse_blastn_check(MAIN_SPECIES, query_orf, sseq_fsa):
                append_to_file(BLASTN_RESULTS_DIR + query_orf, sseq_fsa)
                break
    print("Finished adding " + DB_NAME + " matches to pre-MUSCLE FASTA files")
'''
# RUN MUSCLE
failed_files = []
for multi_fasta_file in glob.glob(os.path.join(BLASTN_RESULTS_DIR + "*")):
    try:
        muscle(multi_fasta_file)
    except:
        failed_files.append(multi_fasta_file)

print('ALL DONE')
print('\n\nTHE FOLLOWING PRE-MUSCLE FILES FAILED')
for i in failed_files:
    print(i)
