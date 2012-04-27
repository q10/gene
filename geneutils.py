import os, glob, subprocess, csv, shelve, datetime
from Bio import SeqIO

DB_DIR = 'db/'
E_VALUE_THRESHOLD = 1e-10
TIME_STACK = []

LOCAL_CDS_DATABASE = dict()
LOCAL_PEP_DATABASE = dict()

def start_timing():
    TIME_STACK.append(datetime.datetime.now())

def finish_timing():
    return datetime.datetime.now() - TIME_STACK.pop()
    
def execute(*commands):
    for comm in commands:
        print(comm)
        if subprocess.call(comm, shell=True, executable="/bin/bash") != 0:
            raise Exception("Command \"" + comm + "\" failed to run!")

def qseq_sseq_pairs(filename):
    for line in csv.reader(open(filename, 'r')):
        yield (line[0], line[1])

def qseq_sseq_sets(filename):
    qs_sets = dict()
    for line in csv.reader(open(filename, 'r')):
        if qs_sets.has_key(line[0]):
            qs_sets[line[0]].append(line[1])
        else:
            qs_sets[line[0]] = [line[1]]
    for qss in qs_sets.iteritems():
        yield qss

def fasta_entries(filename):
    for entry in SeqIO.parse(open(filename),'fasta'):
        yield entry
        
def write_file(filename, contents):
    open(filename, 'w').writelines(contents)

def append_to_file(filename, contents):
    open(filename,'a').writelines(contents)

# generates the BLASTP and BLASTN index files for BLAST to use
def generate_blast_database():
    execute("cd " + DB_DIR + " && rm -rf *.pin *.psq *.phr *.nin *.nsq *.nhr")
    for filename in os.listdir(DB_DIR):
        base_command = "cd " + DB_DIR + " && makeblastdb -in " + filename + " -out " + filename[:-8]
        if filename[-7:] == 'pep.fsa':
            execute(base_command+" -dbtype prot")
        elif filename[-7:] == 'cds.fsa':
            execute(base_command+" -dbtype nucl")

def generate_python_databases(SUBJECT_DBS):
	for DB_NAME in SUBJECT_DBS:
		for entry in fasta_entries(DB_DIR + DB_NAME + '_cds.fsa'):
			LOCAL_CDS_DATABASE[DB_NAME, entry.id] = entry
		for entry in fasta_entries(DB_DIR + DB_NAME + '_pep.fsa'):
			LOCAL_PEP_DATABASE[DB_NAME, entry.id] = entry
		print("added " + DB_NAME + " to PEP and CDS databases")
	
	print("saving CDS and PEP DBs to file...")
	shelf = shelve.open('FAST_DATABASE')
	shelf['LOCAL_CDS_DATABASE'] = LOCAL_CDS_DATABASE
	shelf['LOCAL_PEP_DATABASE'] = LOCAL_PEP_DATABASE
	shelf.close()
	print("done.")

def load_python_databases():
    shelf = shelve.open('FAST_DATABASE')
    LOCAL_CDS_DATABASE = shelf['LOCAL_CDS_DATABASE']
    LOCAL_PEP_DATABASE = shelf['LOCAL_PEP_DATABASE']
    shelf.close()
			
def blastp(query_fasta_filepath, db_namepath, evalue=E_VALUE_THRESHOLD, outfmt="\"10 qseqid sseqid evalue\"", outname=""):
    if outname is "":
        outname = query_fasta_filepath.split('/')[-1] + "-" + db_namepath.split('/')[-1] + ".blastp.csv"
    comm = "blastp -query " + query_fasta_filepath + \
        " -outfmt " + outfmt + \
        " -evalue " + str(evalue) + \
        " -db " + db_namepath + \
        " -out " + outname
    execute(comm)
    print("finished BLASTP of query " + query_fasta_filepath + " against database " + db_namepath)
    
def blastn(query_fasta_filepath, db_namepath, evalue=E_VALUE_THRESHOLD, outfmt="\"10 qseqid sseqid evalue\"", outname=""):
    if outname is "":
        outname = query_fasta_filepath.split('/')[-1] + "-" + db_namepath.split('/')[-1] + ".blastn.csv"
    comm = "blastn -query " + query_fasta_filepath + \
        " -outfmt " + outfmt + \
        " -evalue " + str(evalue) + \
        " -db " + db_namepath + \
        " -out " + outname
    execute(comm)
    print("finished BLASTN of query " + query_fasta_filepath + " against database " + db_namepath)

def muscle(infile, outfile="", clw=False):
    if outfile is "":
        if clw:
            outfile = infile + ".afa"
        else:
            outfile = infile + ".fasta"
    comm = "muscle" + \
        " -in " + infile + \
        " -out " + outfile
    if clw:
        comm += "-clw" # for ClustalW format output instead of FASTA
    execute(comm)
    print("finished running MUSCLE on " + infile)
    
# write seq to short FASTA file, BLAST short fasta file, check results, return boolean
def reverse_blastp_check(orig_qdb, orig_qorf, orig_sseq_fsa):
    write_file('temp_queryp.fasta', orig_sseq_fsa)
    blastp('temp_queryp.fasta', DB_DIR+orig_qdb, outfmt="\"10 sseqid evalue\"", outname="temp_blastp.csv")
    for line in csv.reader(open('temp_blastp.csv', 'r')):
        if line[0] == orig_qorf:
            return True
    return False

def reverse_blastn_check(orig_qdb, orig_qorf, orig_sseq_fsa):
    write_file('temp_queryn.fasta', orig_sseq_fsa)
    blastn('temp_queryn.fasta', DB_DIR+orig_qdb, outfmt="\"10 sseqid evalue\"", outname="temp_blastn.csv")
    for line in csv.reader(open('temp_blastn.csv', 'r')):
        if line[0] == orig_qorf:
            return True
    return False
