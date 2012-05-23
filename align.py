from geneutils import *

DB_DIR = '../Proteomes/Candida_proteomes/'
BLASTN_RESULTS_DIR = 'resultsn/'

#MAIN_SPECIES = 'S288C'
#SUBJECT_DBS = ['AWRI1631_ABSV01000000', 'AWRI796_ADVS01000000', 'CBS7960_AEWL01000000', \
#'CLIB215_AEWP01000000', 'CLIB324_AEWM01000000', 'CLIB382_AFDG01000000', 'EC1118_PRJEA37863', \
#'EC9-8_AGSJ01000000', 'FL100_AEWO01000000', 'FostersB_AEHH01000000', 'FostersO_AEEZ01000000', \
#'Kyokai7_BABQ01000000', 'LalvinQA23_ADVV01000000', 'M22_ABPC01000000', 'PW5_AFDC01000000', \
#'RM11-1a_AAEG01000000', 'Sigma1278b_ACVY01000000', 'JAY291_ACFL01000000', 'T7_AFDE01000000', \
#'T73_AFDF01000000', 'UC5_AFDD01000000', 'Vin13_ADXC01000000', 'VL3_AEJS01000000', 'Y10_AEWK01000000', \
#'YJM269_AEWN01000000', 'YJM789_AAFW02000000', 'YPS163_ABPD01000000', 'W303_MPG_2011']

#SUBJECT_DBS = ['Ashbya_gossypii', 'Candida_glabrata', 'K_waltii', 'Kluyveromyces_lactis', \
#'Kluyveromyces_thermotolerans', 'S_castellii_040406', 'Saccharomyces_kluyveri', 'S_bayanus_ORFs.fasta', \
#'S_mikatae_ORFs.fasta', 'S_paradoxus_ORFs.fasta']

#SUBJECT_DBS = ['Ashbya_gossypii', 'Candida_glabrata', 'K_waltii', 'Kluyveromyces_lactis', \
#'Pichia_sorbitophila', 'S_castellii_040406', 'Saccharomyces_kluyveri', 'Yarrowia_lipolytica', \
#'Zygosaccharomyces_rouxii', 'aspergillus_fumigatus_1', 'aspergillus_nidulans_fgsc_a4_1', \
#'aspergillus_niger_1', 'aspergillus_oryzae_1', 'aspergillus_terreus_1', 'candida_albicans_sc5314_assembly_21_1', \
#'candida_albicans_wo-1_1', 'candida_guilliermondii_1', 'candida_lusitaniae_1', 'candida_parapsilosis_1', \
#'candida_tropicalis_3', 'chaetomium_globosum_1', 'debaryomyces_hansenii_1', 'histoplasma_capsulatum_g186ar_2', \
#'histoplasma_capsulatum_h143_2', 'histoplasma_capsulatum_h88_2', 'histoplasma_capsulatum_nam1_1', \
#'lodderomyces_elongisporus_1', 'magnaporthe_grisea_70-15_mitochondria_6', 'magnaporthe_grisea__m._oryzae__70-15_6', \
#'neosartorya_fischeri', 'neurospora_crassa_or74a__finished__10', 'schizosaccharomyces_cryophilus_oy26_3', \
#'schizosaccharomyces_japonicus_yfs275_4', 'schizosaccharomyces_octosporus_5', 'schizosaccharomyces_pombe_972h-_2', \
#'sclerotinia_sclerotiorum_2', 'uncinocarpus_reesii_2', 'S_bayanus', 'S_mikatae', 'S_paradoxus']


MAIN_SPECIES_LIST = ['S288C_cerevisiae__shift1', 'S288C_cerevisiae__shift2']
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
generate_blast_database()


execute('mkdir -p ' + BLASTN_RESULTS_DIR + ' && rm -rf ' + BLASTN_RESULTS_DIR + '*')
for entry in fasta_entries(DB_DIR + MAIN_SPECIES + '_cds.fsa'):
    write_file(BLASTN_RESULTS_DIR + entry.id, entry.format('fasta'))
    print("created pre-MUSCLE fasta file for " + entry.id)

# RUN FULL BLASTN BETWEEN S288C AND EACH OF THE OTHER STRAINS
for MAIN_SPECIES in MAIN_SPECIES_LIST:
    for DB_NAME in SUBJECT_DBS:
        blast('P', DB_DIR + MAIN_SPECIES + "_pep.fsa", DB_DIR + DB_NAME, evalue=1e-5, outfmt=0, outname=MAIN_SPECIES + "--" + DB_NAME + ".blastp.csv")

print('ALL DONE')