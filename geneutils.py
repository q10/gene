import os, glob, csv, shelve, datetime, subprocess as sp
from Bio import SeqIO


"""
These values are the defaults that will be used.
"""
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
    BASH interpreter.  Also captures stderr output.  Throws exceptions if command executes with error.
    """
    outputs = []
    for comm in commands:
        print(comm)
        outputs.append(sp.check_output(comm, stderr=sp.STDOUT, shell=True, executable="/bin/bash"))
        print(outputs[-1])
    outputs = outputs[0] if len(commands) is 1 else outputs
    return outputs

def detailed_execute(*commands):
    """
    Generic handler that passes a command string or list of command strings into the
    BASH interpreter.  Captures stderr in a separate stream.  Throws exceptions if command executes with error.
    Returns a pair of streams, (<stdout(s)>, <stderr(s)>)
    """
    outs, errors = [], []
    for comm in commands:
        p = sp.Popen(comm, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, executable="/bin/bash")
        stdout, stderr = p.communicate()
        outs.append(stdout), errors.append(errors)
        print(outs[-1])
    outs = outs[0] if len(commands) is 1 else outs
    errors = errors[0] if len(commands) is 1 else errors
    return (outs, errors)

def list2string(*strvars):
    """
    Converts a list of arguments of any type into string form.  Written for "easier" concatenation of variables and strings
    """
    return 'ls'
    #a = "".join(map(str, strvars))
    #print('THIS IS A')
    #print(a)
    #return a

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

def generate_blast_database(DB_DIR):
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

def generate_python_pep_database(DB_DIR, SUBJECT_DBS):
    """
    Loads all entries in specified BLASTP databases into memory as FASTA objects.  This is for Python's use, and not for BLAST*'s use.
        - *SUBJECT_DBS: list of database names (without filename extensions).  The function will automatically look for the database in DB_DIR directory
        Fasta files must be named "<DB-NAME>_pep.fsa"
    """
    for DB_NAME in SUBJECT_DBS:
        for entry in fasta_entries(DB_DIR + DB_NAME + '_pep.fsa'):
            LOCAL_PEP_DATABASE[DB_NAME, entry.id] = entry
        print("added " + DB_NAME + " to PEP database")


def generate_python_cds_database(DB_DIR, SUBJECT_DBS):
    """
    Loads all entries in specified BLASTP databases into memory as FASTA objects.  This is for Python's use, and not for BLAST*'s use.
        - *SUBJECT_DBS: list of database names (without filename extensions).  The function will automatically look for the database in DB_DIR directory
        Fasta files must be named "<DB-NAME>_cds.fsa"
    """
    for DB_NAME in SUBJECT_DBS:
        for entry in fasta_entries(DB_DIR + DB_NAME + '_cds.fsa'):
            LOCAL_CDS_DATABASE[DB_NAME, entry.id] = entry
        print("added " + DB_NAME + " to CDS database")
	
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

def blast(btype, query_fasta_filepath, db_namepath, evalue=E_VALUE_THRESHOLD, outfmt="\"10 qseqid sseqid evalue\"", outname=""):
    """
    Runs a BLAST* query.
        - btype: 'N' or 'P' - specifies nucleotide or protein sequence
        - query_fasta_filepath: the file containing the query entry in FASTA format
        - db_namepath: the path to database + name of the database (WITHOUT the filename extensions)
        - evalue=E_VALUE_THRESHOLD: the maximum expectation value threshold to use
        - outfmt="\"10 qseqid sseqid evalue\"": the format of the query output.  Default is a 3-column file, containing the query ID, subject ID, and e-value
        - outname="": the name of the output file.  By default, it will be "<QUERY-FILENAME>-<SUBJECT-DB-NAME>.blast<btype>.csv"
    """
    exe, ext = ("blastn", ".blastn.csv") if btype is 'N' or btype is 'n' else ("blastp", ".blastp.csv")
    if outname is "":
        outname = list2string(query_fasta_filepath.split('/')[-1], '--', db_namepath.split('/')[-1], ext)
    execute(list2string(exe, " -query ", query_fasta_filepath, " -outfmt ", outfmt, " -evalue ", evalue, " -db " + db_namepath, " -out ", outname))
    #print(list2string("finished BLAST", btype, " of query ", query_fasta_filepath, " against database ", db_namepath))

def muscle(infile, outfile="", clw=False):
    """
    Runs a MUSCLE query.
        - infile: the name of the input file.
        - outfile="": the name of the output file.  Default will be "<infile>.fasta"
        - clw=False: flag for ClustalW format output instead of FASTA
    """
    if outfile is "":
        outfile = infile + ".afa" if clw else infile + ".fasta"
    comm = list2string("muscle", " -in ", infile, " -out ", outfile)
    comm = comm + "-clw" if clw else comm
    execute(comm)
    print("finished running MUSCLE on " + infile)

# write seq to short FASTA file, BLAST short fasta file, check results, return boolean
def reverse_blast_check(btype, orig_qdb, orig_qorf, orig_sseq_fsa, evalue=E_VALUE_THRESHOLD):
    """
    Performs a reverse BLAST* check to assert that a BLAST* hit is valid
        - btype: 'N' or 'P' - specifies nucleotide or protein sequence
        - orig_qdb: the ORF name of the candidate subject sequence
        - orig_qorf: the ORF name of the query sequence
        - orig_sseq_fsa: the candidate subject sequence (in FASTA format) that is to be reverse-checked against the original query sequence 
        - evalue=E_VALUE_THRESHOLD: the maximum expectation value threshold to use
    """
    tmpq, outf = ('temp_queryn.fasta', 'temp_blastn.csv') if btype is 'N' or btype is 'n' else ('temp_queryp.fasta', 'temp_blastp.csv')
    write_file(tmpq, orig_sseq_fsa)
    blast(btype, tmpq, orig_qdb, evalue, outfmt="\"10 sseqid evalue\"", outname=outf)
    return any(csv.reader(open(outf, 'r')), lambda line: line[0] == orig_qorf)

def frame_shift(fas, shiftamt):
    """
    Returns a FASTA entry that is frame-shifted by the specified amount (1 or 2 letters)
    Renames the FASTA ID to reflect the change
        - fas: the FASTA entry object
        - shiftamt: the shift amount.  Must be either 1 or 2, or exception will be thrown
    """
    if shiftamt != 1 and shiftamt != 2:
        raise Exception('FASTA frame-shifting must be by 1 or 2')
    fas.id = fas.id + '--SHIFT_' + str(shiftamt)
    fas.name = fas.name + '--SHIFT_' + str(shiftamt)
    fas.seq = fas.seq.lstrip(fas.seq[0:shiftamt])
    return fas

def fasta_translate(fas):
    """
    Returns a FASTA entry that is translated into the corresponding protein sequence
    Renames the FASTA ID to reflect the change
        - fas: the FASTA entry object.  Must contain nucleotide sequences
    """
    fas.seq = fas.seq.translate()
    return fas

def generate_fasta_frame_shifted_file(filename, shiftamt):
    """
    Creates a new FASTA file where all the entries from the original FASTA file are 
    frame-shifted by the amount specified.  Each FASTA entry will be renamed to reflect this change.
    The new file's name will reflect this change and the old file remains unchanged.
        - filename: the name of the FASTA file to be frame-shifted
        - shiftamt: the shift amount.  Must be either 1 or 2, or exception will be thrown
    """
    if shiftamt != 1 and shiftamt != 2:
        raise Exception('FASTA file frame-shifting must be by 1 or 2')
    outfile = os.path.splitext(filename)[0] + "__shift" + str(shiftamt) + ".fasta"
    for entry in fasta_entries(filename):
        append_to_file(outfile, frame_shift(entry, shiftamt).format('fasta'))
        
def generate_fasta_translated_file(filename):
    """
    Creates a new FASTA file where all the entries from the original FASTA file are 
    translated into the corresponding protein sequences.  Each FASTA entry will be renamed to reflect this change.
    The new file's name will reflect this change and the old file remains unchanged.
        - filename: the name of the FASTA file to be frame-shifted.  Must contain nucleotide sequences
    """
    outfile = os.path.splitext(filename)[0] + "__translated" + ".fasta"
    for entry in fasta_entries(filename):
        append_to_file(outfile, fasta_translate(entry, shiftamt).format('fasta'))

def generate_fasta_frame_shifted_translated_file(filename, shiftamt):
    """
    Creates a new FASTA file where all the entries from the original FASTA file are 
    frame-shifted by the amount specified, and then translated to protein sequence.
    Each FASTA entry will be renamed to reflect this change.
    The new file's name will reflect this change and the old file remains unchanged.
        - filename: the name of the FASTA file to be frame-shifted.  Must contain nucleotide sequences
        - shiftamt: the shift amount.  Must be either 1 or 2, or exception will be thrown
    """
    if shiftamt != 1 and shiftamt != 2:
        raise Exception('FASTA file frame-shifting must be by 1 or 2')
    outfile = os.path.splitext(filename)[0] + "__shift" + str(shiftamt) + "__translated" + ".fasta"
    for entry in fasta_entries(filename):
        append_to_file(outfile, fasta_translate(frame_shift(entry, shiftamt)).format('fasta'))
