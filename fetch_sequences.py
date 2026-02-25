"""
fetch_sequences.py ------------------------------------------------------

Author: Paul Antezana

This serves to create a python script to individually fetch 
and download sequences for both the Bowhead and Blue whales. 

It then saves them, as well as data about their % identity 
into the raw folder, which will be used for analysis later.

-------------------------------------------------------------------------

Gene sequences for the bowhead and blue whale:
 ---> Blue whale: downloaded from NCBI (GenBank)
 ---> Bowhead whale: extracted from the Bowhead Whale Genome Resource

### BEFORE RUNNING:##################
  Bowhead CDS file must be 
  downloaded from bowhead-whale.org.
#####################################
  
"""

# import all the modules necessary for this program
# these were found to be abundant in many real bioinformatics projects
import os
import time
import subprocess
import shutil

# Entrez specifically is important because it allows us to implement code that helps us directly access the NCBI database
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
from io import StringIO

# NCBI requires an email address to use their API
Entrez.email = "INSERT EMAIL ADDRESS HERE" # <----------------------------------------- IMPORTANT! Insert your email here first.

# list the 5 genes we are comparing between species
geneList = ["CIRBP", "RPA2", "PCNA", "ERCC1", "TP53"]

# build the output folder path relative to this script, then create it (if needed)
script_dir = os.path.dirname(os.path.abspath(__file__))

# now place the raw data into this path
RAW_FOLDER = os.path.join(script_dir, "../data/raw")
os.makedirs(RAW_FOLDER, exist_ok=True)

# this variable stores the path information for the gene file download from the Bowhead Whale Genome Resource website.
# it's important to know that manually downloading this fasta file is important because the bowhead whale genome
# cannot be accessed from GenBank as it's unique to only this website.
BOWHEAD_CDS = "C:/Users/paula/Downloads/bowhead_whale_coding_sequences/bowhead_whale_coding_sequences.fasta"

# this varaible stores the path information for the BLAST path
BLAST_BIN = "C:/Users/paula/OneDrive/UWB/Classes/Winter Quarter 2026/Bioinformatics/blast-2.17.0+/bin"

# this temporary folder is made to create no spaces in the path. BLAST+ seems to have an issue with spaces
# in paths and as for now, this just serves as a temporary solution until a better one is written.
TEMP_FOLDER = "C:/Users/paula/Downloads/bowhead_blast_temp/"

# minimum identity % to accept a BLAST hit as a real match. 95% might be too high to hit anything, so 70% is a good lowball
# estimate.
MIN_IDENTITY = 70.0

def find_gene_id(gene_name):
    """
    Purpose: Search NCBI's Gene database for a gene by name in the blue whale genome.

    Parameters:
        gene_name (str): the gene symbol to search for (e.g. "CIRBP")

    Returns:
        str: the NCBI Gene ID if found, or None if no match exists
    """
    query = f"{gene_name}[Gene Name] AND txid9771[Organism:exp]"
    handle = Entrez.esearch(db="gene", term=query, retmax=5)
    record = Entrez.read(handle)
    handle.close()
    # if it found something, return the first matching gene ID. the database automatically lists then
    # from ascending order, so the top one will have the highest % identity.
    if record["IdList"]:
        return record["IdList"][0]
    # if it didn't, return nothing.
    return None

def get_mrna_accession(gene_id):
    """
    Purpose: Fetch the gene record from NCBI and extract the mRNA accession number.

    Parameters:
        gene_id (str): the NCBI Gene ID returned by find_gene_id()

    Returns:
        str: the mRNA accession number (e.g. "NM_001234") if found
        None: if mRNA accession number is NOT found
    """
    handle = Entrez.efetch(db="gene", id=gene_id, rettype="gene_table", retmode="text")
    text = handle.read()
    handle.close()
    for line in text.split("\n"):
        for word in line.split():
            # mRNA accession usually starts with either NM or XM_
            if word.startswith("NM_") or word.startswith("XM_"):
                # remove any trailing punctuation
                return word.strip(".,;")
    return None

def fetch_fasta(accession):
    """
    Purpose: Download the DNA sequence for a given accession number from NCBI in FASTA format.

    Parameters:
        accession (str): the mRNA accession number to download (e.g. "NM_XXXXX")

    Returns:
        str: the full FASTA-formatted sequence as a text string
    """
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    sequence = handle.read()
    handle.close()
    return sequence

print("### Step 1: Fetching blue whale sequences from NCBI ###\n")

# for EACH gene in the geneList (storing our genes of interest), find the gene ID, get the mRNA
# accession number from it, and then download the FASTA sequence. It then saves it.
for gene in geneList:
    print(f"  {gene}...", end=" ")

    # loop up the given gene ID
    gene_id = find_gene_id(gene)
    accession = None

    # if there IS a gene ID, return the mRNA accession number
    if gene_id:
        accession = get_mrna_accession(gene_id)
    # otherwise, print out an error message.
    if accession is None:
        print("not found in NCBI.")
        continue

    # download the FASTA file and put it into our RAW_FOLDER file.
    fasta_text = fetch_fasta(accession)
    out_file = os.path.join(RAW_FOLDER, f"{gene}_blue_whale.fasta")
    with open(out_file, "w") as f:
        f.write(fasta_text)

    # After the entire code runs for a given gene, wait 1 second. This is to comply with NCBI's usage policy
    # otherwise it may not work.
    print(f"saved ({accession})")
    time.sleep(1)

print("\n### Step 2: Extracting bowhead sequences from bowhead-whale.org CDS file ###\n")

'''
    There is a problem here. Since the Bowhead Whale genome can only be accessed in the bowhead
    whale genome resource website, we cannot use our previous code to access it. We must use the local
    BLAST we have and find the bowhead version of each gene within geneList by comparing them against
    the blue whale sequence we just recently downloaded.
'''

os.makedirs(TEMP_FOLDER, exist_ok=True)

# this part just adds the BLAST+ bin to PATH so we can start using makeblastdb and blastn
# calls by name alone. Otherwise, it can be tedious to call them.
env = os.environ.copy()
env["PATH"] = BLAST_BIN + ";" + env["PATH"]

db_path = TEMP_FOLDER + "bowhead_cds_db"

# build the BLAST database from the bowhead CDS file (only needs to be done once per run of these programs)

# this checks if it already exists. if it does, run the following code.
if not os.path.exists(db_path + ".nsq") and not os.path.exists(db_path + ".nsi"):
    print("Building local BLAST database from bowhead file...")
    cmd = ["makeblastdb", "-in", BOWHEAD_CDS, "-dbtype", "nucl",
           "-out", db_path, "-parse_seqids"]
    
    # store the results into the result variable
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)

    # if the return code is not equal to 0, we print out an error message
    if result.returncode != 0:
        print(f"ERROR building database:\n{result.stderr}")
        exit(1)
    print("Database built.\n")
else:
    # it does exist and therefore print out a helpful message
    print("BLAST database already exists — skipping build.\n")

# now for each gene in the geneList, this will copy over the blue whale file to the temporary folder with NO spaces.
# this is a fix in response to the problem BLAST seems to have if it does have spaces in paths.
for gene in geneList:
    blue_src = os.path.join(RAW_FOLDER, f"{gene}_blue_whale.fasta")
    if not os.path.exists(blue_src):
        print(f"{gene:<10} SKIPPED — no blue whale sequence to use as query")
        continue

    blue_tmp = TEMP_FOLDER + f"{gene}_query.fasta"
    blast_out = TEMP_FOLDER + f"{gene}_blast.xml"
    shutil.copy2(blue_src, blue_tmp)

    print(f"{gene:<10} BLAST-ing...", end=" ", flush=True)

    # this part runs blastn with the query set as the blue whale sequence and the database set as the bowhead 
    # coding DNA sequence.

    # also note, we're setting the e-value to be very small so the results are filtered in a way that
    # only statistically significant results show up.
    cmd = ["blastn", "-query", blue_tmp, "-db", db_path, "-out", blast_out,
           "-outfmt", "5", "-max_target_seqs", "5", "-evalue", "1e-10", "-num_threads", "4"]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)

    if result.returncode != 0:
        print(f"BLAST error: {result.stderr.strip()}")
        continue

    # read the BLAST output file
    with open(blast_out) as f:
        xml_text = f.read()

    # if the BLAST output is missing, print out an error message
    if not xml_text.strip() or "<BlastOutput" not in xml_text:
        print("no BLAST output.")
        continue

    # parse the XML results
    blast_record = NCBIXML.read(StringIO(xml_text))

    # if there are no known alignments in the blast database, then print out an error message
    if not blast_record.alignments:
        print("no hits found in bowhead CDS.")
        continue

    # BLAST automatically lists the most accurate % identity hit at the very start of the list, so store that first result
    # into top_align.
    top_align = blast_record.alignments[0]
    top_hsp = top_align.hsps[0]

    # calculate percent identity and round to 2 decimal places
    identity = round((top_hsp.identities / top_hsp.align_length) * 100, 2)
    # internal bowhead sequence ID (e.g. bmy_XXXXX)
    hit_id = top_align.hit_id  
    evalue = top_hsp.expect

    # if the given hit is below the MINUMUM identity (i.e. is too different to be a real ortholog), then print out the identity
    # and minimum identity in a message.
    if identity < MIN_IDENTITY:
        print(f"top hit only {identity}% identity (below {MIN_IDENTITY}%) — skipping.")
        continue

    # finds and retrieves the full bowhead sequence from the CDS file using its ID.
    bowhead_seq = None
    for record in SeqIO.parse(BOWHEAD_CDS, "fasta"):
        # if the record id matches the current hit id, then save the sequence into the bowhead_seq variable
        if record.id == hit_id:
            bowhead_seq = record
            break
    # if the bowhead sequence could not be retrieved, then print out a message
    if bowhead_seq is None:
        print(f"could not retrieve sequence for {hit_id}.")
        continue

    # save the bowhead sequence to the raw folder (in data so dir should be data/raw)
    out_file = os.path.join(RAW_FOLDER, f"{gene}_bowhead.fasta")
    with open(out_file, "w") as f:
        f.write(f">{gene}_bowhead_whale [source:{hit_id}] [identity:{identity}%] [evalue:{evalue:.2e}]\n")
        f.write(str(bowhead_seq.seq) + "\n")

    # print out the gene's % identity and its given hit ID
    print(f"saved! ({identity}% identity, source: {hit_id})")

# if the program ran without no major errors, print out this message. this serves to verify everything went ok.
print("\nDone! Sequences saved to data/raw/")
print("###################################################")
print("NEXT STEPS: Run blast_compare to compare results.")
print("###################################################")
