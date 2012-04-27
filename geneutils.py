import os, glob, subprocess, csv, shelve, datetime
from Bio import SeqIO


"""
These values must be set prior to using the library, or defaults will be used.
"""
DB_DIR = 'db/'
E_VALUE_THRESHOLD = 1e-10


"""
These variables are internal to the library.
"""
TIME_STACK = []
LOCAL_CDS_DATABASE = dict()
LOCAL_PEP_DATABASE = dict()

def start_timing():
    """
    Starts timing execution from evocation point.  May be called sequentially for sub-timing.
    """
    TIME_STACK.append(datetime.datetime.now())

def finish_timing():
    """
    Finishes timing execution since the last start_timing() call, and
    returns the time difference.
    """
    return datetime.datetime.now() - TIME_STACK.pop()
    
def execute(*commands):
    """
    Generic handler that passes a command string or list of command strings into the
    BASH interpreter.  Throws exceptions if command executes with error.
    """
    for comm in commands:
        print(comm)
        if subprocess.call(comm, shell=True, executable="/bin/bash") != 0:
            raise Exception("Command \"" + comm + "\" failed to run!")

def qseq_sseq_pairs(filename):
    """
    Iterates over a text file and yields pairs of values, where each pair corresponds to the values in the 
    first and second column of each line.
        - filename: the full name of the file.  Must be CSV-formatted with at least 2 columns of text
    """
    for line in csv.reader(open(filename, 'r')):
        yield (line[0], line[1])

def qseq_sseq_sets(filename):
    """
    Iterates over a BLAST* results file with the format "\"10 qseqid sseqid evalue\"" and 
    yields dictionaries, where each dictionary contains all the candidate subject ORFS, ordered by 
    increasing e-value, that match one query ORF in the FASTA query file.
        - filename: the full name of the file.  Must be CSV-formatted with at least 2 columns of text
    """
    qs_sets = dict()
    for line in csv.reader(open(filename, 'r')):
        if qs_sets.has_key(line[0]):
            qs_sets[line[0]].append(line[1])
        else:
            qs_sets[line[0]] = [line[1]]
    for qss in qs_sets.iteritems():
        yield qss

def fasta_entries(filename):
    """
    Iterates over all entries in a standard FASTA file. 
    """
    for entry in SeqIO.parse(open(filename),'fasta'):
        yield entry
        
def write_file(filename, contents):
    """
    Writes contents to file.  Creates a new file if nonexistent, and overwrites existent files.
    """
    open(filename, 'w').writelines(contents)

def append_to_file(filename, contents):
    """
    Appends contents to existent file.  Creates a new file if nonexistent.
    """
    open(filename,'a').writelines(contents)

def generate_blast_database():
    """
    Goes into the directory specified by DB_DIR, and generates the BLAST* index files from the FASTA files for BLAST* to use.
    As a precaution, removes all old BLAST* index files for clean rebuild of BLAST databases.
    Need to call before main program runs, or make sure BLAST* index files are in place prior to running.
    """
    execute("cd " + DB_DIR + " && rm -rf *.pin *.psq *.phr *.nin *.nsq *.nhr")
    for filename in os.listdir(DB_DIR):
        base_command = "cd " + DB_DIR + " && makeblastdb -in " + filename + " -out " + filename[:-8]
        if filename[-7:] == 'pep.fsa':
            execute(base_command+" -dbtype prot")
        elif filename[-7:] == 'cds.fsa':
            execute(base_command+" -dbtype nucl")

def generate_python_databases(*SUBJECT_DBS):
    """
    Loads all entries in specified BLASTP databases into memory as FASTA objects.  This is for Python's use, and not for BLAST*'s use.
        - *SUBJECT_DBS: list of database names (without filename extensions).  The function will automatically look for the database in DB_DIR directory 
    """
	for DB_NAME in SUBJECT_DBS:
		for entry in fasta_entries(DB_DIR + DB_NAME + '_cds.fsa'):
			LOCAL_CDS_DATABASE[DB_NAME, entry.id] = entry
		for entry in fasta_entries(DB_DIR + DB_NAME + '_pep.fsa'):
			LOCAL_PEP_DATABASE[DB_NAME, entry.id] = entry
		print("added " + DB_NAME + " to PEP and CDS databases")
	
	#print("saving CDS and PEP DBs to file...")
	#shelf = shelve.open('FAST_DATABASE')
	#shelf['LOCAL_CDS_DATABASE'] = LOCAL_CDS_DATABASE
	#shelf['LOCAL_PEP_DATABASE'] = LOCAL_PEP_DATABASE
	#shelf.close()
	#print("done.")

#def load_python_databases():
#    shelf = shelve.open('FAST_DATABASE')
#    LOCAL_CDS_DATABASE = shelf['LOCAL_CDS_DATABASE']
#    LOCAL_PEP_DATABASE = shelf['LOCAL_PEP_DATABASE']
#    shelf.close()
			
def blastp(query_fasta_filepath, db_namepath, evalue=E_VALUE_THRESHOLD, outfmt="\"10 qseqid sseqid evalue\"", outname=""):
    """
    Runs a BLASTP query.
        - query_fasta_filepath: the file containing the query entry in FASTA format
        - db_namepath: the path to database + name of the database (WITHOUT the filename extensions)
        - evalue=E_VALUE_THRESHOLD: the maximum expectation value threshold to use
        - outfmt="\"10 qseqid sseqid evalue\"": the format of the query output.  Default is a 3-column file, containing the query ID, subject ID, and e-value
        - outname="": the name of the output file.  By default, it will be "<QUERY-FILENAME>-<SUBJECT-DB-NAME>.blastp.csv"
    """
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
    """
    Runs a BLASTN query.
        - query_fasta_filepath: the file containing the query entry in FASTA format
        - db_namepath: the path to database + name of the database (WITHOUT the filename extensions)
        - evalue=E_VALUE_THRESHOLD: the maximum expectation value threshold to use
        - outfmt="\"10 qseqid sseqid evalue\"": the format of the query output.  Default is a 3-column file, containing the query ID, subject ID, and e-value
        - outname="": the name of the output file.  By default, it will be "<QUERY-FILENAME>-<SUBJECT-DB-NAME>.blastn.csv"
    """
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
    """
    Runs a MUSCLE query.
        - infile: the name of the input file.
        - outfile="": the name of the output file.  Default will be "<infile>.fasta"
        - clw=False: flag for ClustalW format output instead of FASTA
    """
    if outfile is "":
        if clw:
            outfile = infile + ".afa"
        else:
            outfile = infile + ".fasta"
    comm = "muscle" + \
        " -in " + infile + \
        " -out " + outfile
    if clw:
        comm += "-clw"
    execute(comm)
    print("finished running MUSCLE on " + infile)
    
# write seq to short FASTA file, BLAST short fasta file, check results, return boolean
def reverse_blastp_check(orig_qdb, orig_qorf, orig_sseq_fsa):
    """
    Performs a reverse BLASTP check to assert that a BLASTP hit is valid
        - orig_qdb: the ORF name of the candidate subject sequence
        - orig_qorf: the ORF name of the query sequence
        - orig_sseq_fsa: the candidate subject sequence (in FASTA format) that is to be reverse-checked against the original query sequence 
    """
    write_file('temp_queryp.fasta', orig_sseq_fsa)
    blastp('temp_queryp.fasta', DB_DIR+orig_qdb, outfmt="\"10 sseqid evalue\"", outname="temp_blastp.csv")
    for line in csv.reader(open('temp_blastp.csv', 'r')):
        if line[0] == orig_qorf:
            return True
    return False

def reverse_blastn_check(orig_qdb, orig_qorf, orig_sseq_fsa):
    """
    Performs a reverse BLASTN check to assert that a BLASTP hit is valid
        - orig_qdb: the ORF name of the candidate subject sequence
        - orig_qorf: the ORF name of the query sequence
        - orig_sseq_fsa: the candidate subject sequence (in FASTA format) that is to be reverse-checked against the original query sequence 
    """
    write_file('temp_queryn.fasta', orig_sseq_fsa)
    blastn('temp_queryn.fasta', DB_DIR+orig_qdb, outfmt="\"10 sseqid evalue\"", outname="temp_blastn.csv")
    for line in csv.reader(open('temp_blastn.csv', 'r')):
        if line[0] == orig_qorf:
            return True
    return False
