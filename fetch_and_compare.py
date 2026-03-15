"""
fetch_and_compare.py -----------------------------------------------

Author: Paul Antezana

This script fetches gene sequences for four species and computes pairwise
identity for two species pairs:

Species:
    Cetaceans: Blue whale vs Bowhead whale
    Non-cetaceans (Suidae): Wild boar vs Common warthog

Genes of interest:
    Cancer: CIRBP, RPA2, PCNA, ERCC1, TP53
    Non-cancer (Housekeeping): GAPDH, ACTB, RPL13A, HPRT1, TBP

It compares both cancer-related genes and housekeeping control genes,
then saves all results to a CSV for statistical analysis.

Species info:
    Blue whale (Balaenoptera musculus)      — NCBI txid9771
    Bowhead whale (Balaena mysticetus)      — bowhead-whale.org CDS file
    Wild boar (Sus scrofa)                  — NCBI txid9823
    Common warthog (Phacochoerus africanus) — NCBI txid41426

### BEFORE RUNNING:##################
  Bowhead CDS file must be
  downloaded from bowhead-whale.org.

  Blast must be downloaded

  Insert proper email for API
  accession

  install necessary modules

  The folder for running this program should look like this:

  [Whale comparative genome analysis folder]
  -> scripts
        -> fetch_and_compare.py
        -> statistical_comparison.py
  -> blast-2.17.0+
  -> bowhead_whale_coding_sequences.fasta
  
#####################################

"""

# import all the modules necessary for this program
import os
import csv
import time
import subprocess
import shutil
import tempfile

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIXML
from Bio.Align import PairwiseAligner
from io import StringIO

# NCBI requires an email address to use their API
# (Just for purposes of this class, I'll leave my e-mail here)
Entrez.email = "paulag24@uw.edu" # <-------- ### IMPORTANT: Please input a valid email here for the program to work ###

# build file paths relative to this script's location
# important to note that blast freaks out when there are spaces in its path, so 
# the temp_folder is there to fix this issue
script_dir = os.path.dirname(os.path.abspath(__file__))
RAW_FOLDER = os.path.join(script_dir, "../data/raw")
RESULTS_FOLDER = os.path.join(script_dir, "../results")
BOWHEAD_CDS = os.path.join(script_dir, "../bowhead_whale_coding_sequences.fasta")
TEMP_FOLDER = os.path.join(tempfile.gettempdir(), "bowhead_blast_temp")
BLAST_BIN = os.path.join(script_dir, "../blast-2.17.0+/bin")
os.makedirs(RAW_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)
os.makedirs(TEMP_FOLDER, exist_ok=True)

### User can put their genes of interest here and test it without heavily modifying the programs ###

# cancer-related genes
CANCER_GENES = ["CIRBP", "RPA2", "PCNA", "ERCC1", "TP53"]

# housekeeping/control genes
CONTROL_GENES = ["GAPDH", "ACTB", "RPL13A", "HPRT1", "TBP"]

ALL_GENES = CANCER_GENES + CONTROL_GENES

# species taxonomy IDs for NCBI searches
# NOTICE: the bowhead whale is NOT on NCBI so it is handled separately with local BLAST
SPECIES = {
    "blue_whale": "txid9771",
    "wild_boar": "txid9823",
    "warthog": "txid41426",
}

# minimum identity % to accept a BLAST hit as a real match
MIN_IDENTITY = 70.0

### FUNCTIONS ###

def find_gene_id(gene_name, taxon_id):
    """
    Purpose: Search NCBI's gene database for a gene in a specific organism.

    Parameters:
        gene_name (str): gene symbol (e.g. "CIRBP")
        taxon_id (str): NCBI taxonomy ID (e.g. "txid9771" for blue whale)

    Returns:
        str: the NCBI Gene ID if found
        None: if the NCBI Gene ID is not found
    """
    query = f"{gene_name}[Gene Name] AND {taxon_id}[Organism:exp]"
    handle = Entrez.esearch(db="gene", term=query, retmax=5)
    record = Entrez.read(handle)
    handle.close()
    # since NBCI returns the highest result at the top of the list, we search for the result in index 0
    if record["IdList"]:
        return record["IdList"][0]
    return None

def get_mrna_accession(gene_id):
    """
    Purpose: Fetch the gene record and extract the mRNA accession number.

    Parameters:
        gene_id (str): the NCBI Gene ID

    Returns:
        str: mRNA accession (e.g. "NM_001234") or None
    """
    handle = Entrez.efetch(db="gene", id=gene_id, rettype="gene_table", retmode="text")
    text = handle.read()
    handle.close()
    for line in text.split("\n"):
        for word in line.split():
            # typically, mRNA accession usually starts with either NM_ or XM_
            if word.startswith("NM_") or word.startswith("XM_"):
                return word.strip(".,;")
    return None

def fetch_fasta(accession):
    """
    Purpose: Download FASTA sequence from NCBI for a given accession.

    Parameters:
        accession (str): mRNA accession numbr

    Returns:
        str: FASTA-formatted text
    """
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    sequence = handle.read()
    handle.close()
    return sequence

def read_sequence(filepath):
    """
    Purpose: Read a FASTA file and return the DNA sequence as a string.

    Parameters:
        filepath (str): path to the FASTA file

    Returns:
        str: DNA sequence, or None if the file doesn't exist
    """
    if not os.path.exists(filepath):
        return None
    record = SeqIO.read(filepath, "fasta")
    return str(record.seq)

def pairwise_identity(seq1, seq2):
    """
    Purpose: Align two DNA sequences and calculate percent identity.

    Parameters:
        seq1 (str): first DNA sequence
        seq2 (str): second DNA sequence

    Returns:
        float: percent identity, rounded to 2 decimal places
    """
    aligner = PairwiseAligner()
    aligner.mode = "local"
    best = next(iter(aligner.align(seq1, seq2)))
    matches = sum(a == b for a, b in zip(*best))
    return round((matches / best.length) * 100, 2)


### MAIN ###

print("=" * 60)
print("Fetching sequences from NCBI for 3 species")
print("=" * 60)

# fetch blue whale, wild boar, and warthog sequences from NCBI
### REVIEW THIS PART AGAIN AAAAAA ###
for species_label, taxon_id in SPECIES.items():
    new_name = species_label.replace("_", " ").title()
    print(f"\n--- {new_name} ({taxon_id}) ---")

    for gene in ALL_GENES:
        print(f"  {gene}...", end=" ")

        # check if we already downloaded this file (if we did, skip it)
        out_file = os.path.join(RAW_FOLDER, f"{gene}_{species_label}.fasta")
        if os.path.exists(out_file):
            print("already downloaded, skipping.")
            continue

        # look up the gene ID in NCBI for this species
        gene_id = find_gene_id(gene, taxon_id)
        accession = None

        # if we found a gene ID, grab the mRNA accession
        if gene_id:
            accession = get_mrna_accession(gene_id)

        # if we couldn't find an accession, skip this gene
        if accession is None:
            print("not found in NCBI.")
            continue

        # download the FASTA sequence and save it
        fasta_text = fetch_fasta(accession)
        with open(out_file, "w") as f:
            f.write(fasta_text)

        print(f"saved ({accession})")

        # wait 1 second to comply with NCBI usage policy
        time.sleep(1)

### Extract bowhead whale sequences using local BLAST to identify them using blue whale ###

print("\n" + "=" * 60)
print("Extracting bowhead whale sequences via local BLAST")
print("=" * 60)

os.makedirs(TEMP_FOLDER, exist_ok=True)

# add BLAST+ bin to PATH so we can call makeblastdb and blastn directly
env = os.environ.copy()
env["PATH"] = BLAST_BIN + ";" + env["PATH"]
db_path = os.path.join(TEMP_FOLDER, "bowhead_cds_db")

# build BLAST database if it doesn't already exist
if not os.path.exists(db_path + ".nsq") and not os.path.exists(db_path + ".nsi"):
    print("\nBuilding local BLAST database from bowhead CDS file...")
    # copy CDS file to temp folder to avoid spaces-in-path issues with BLAST
    cds_tmp = os.path.join(TEMP_FOLDER, "bowhead_cds.fasta")
    if not os.path.exists(cds_tmp):
        shutil.copy2(BOWHEAD_CDS, cds_tmp)
    cmd = ["makeblastdb", "-in", cds_tmp, "-dbtype", "nucl",
           "-out", db_path, "-parse_seqids"]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        print(f"Error while building database:\n{result.stderr}")
        exit(1)
    print("Database built.\n")
else:
    print("\nBLAST database already exists - skipping build.\n")

# BLAST each blue whale sequence against bowhead CDS to find the bowhead version
for gene in ALL_GENES:
    # use the blue whale sequence as the query
    blue_src = os.path.join(RAW_FOLDER, f"{gene}_blue_whale.fasta")
    if not os.path.exists(blue_src):
        print(f"  {gene:<10} SKIPPED - no blue whale sequence to use as query")
        continue

    # check if bowhead file already exists (skip if so)
    bowhead_out = os.path.join(RAW_FOLDER, f"{gene}_bowhead.fasta")
    if os.path.exists(bowhead_out):
        print(f"  {gene:<10} already extracted, skipping.")
        continue

    # copy query to the temporary folder (avoids spaces-in-path issue with BLAST because it 
    # tends to break if there are spaces in its path)
    blue_tmp = os.path.join(TEMP_FOLDER, f"{gene}_query.fasta")
    blast_out = os.path.join(TEMP_FOLDER, f"{gene}_blast.xml")
    shutil.copy2(blue_src, blue_tmp)

    print(f"  {gene:<10} BLAST-ing...", end=" ", flush=True)

    # run blastn with really strict e-value filters
    cmd = ["blastn", "-query", blue_tmp, "-db", db_path, "-out", blast_out,
           "-outfmt", "5", "-max_target_seqs", "5", "-evalue", "1e-10", "-num_threads", "4"]
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)

    # print error message
    if result.returncode != 0:
        print(f"BLAST error: {result.stderr.strip()}")
        continue

    # read and parse the BLAST output XML
    with open(blast_out) as f:
        xml_text = f.read()
    if not xml_text.strip() or "<BlastOutput" not in xml_text:
        print("no BLAST output.")
        continue
    blast_record = NCBIXML.read(StringIO(xml_text))
    if not blast_record.alignments:
        print("no hits found in bowhead CDS.")
        continue

    # grab the top hit (BLAST lists best matches first)
    top_align = blast_record.alignments[0]
    top_hsp = top_align.hsps[0]
    identity = round((top_hsp.identities / top_hsp.align_length) * 100, 2)
    hit_id = top_align.hit_id
    evalue = top_hsp.expect

    # skip if identity is too low to be a real ortholog based on user-set variable MIN_IDENTITIY
    if identity < MIN_IDENTITY:
        print(f"top hit only {identity}% (below {MIN_IDENTITY}%) — skipping.")
        continue

    # find the full bowhead sequence from the CDS file using its ID
    bowhead_seq = None
    for record in SeqIO.parse(BOWHEAD_CDS, "fasta"):
        if record.id == hit_id:
            bowhead_seq = record
            break
    if bowhead_seq is None:
        print(f"could not retrieve sequence for {hit_id}.")
        continue

    # save the bowhead sequence with appropriate information in the header
    with open(bowhead_out, "w") as f:
        f.write(f">{gene}_bowhead_whale [source:{hit_id}] [identity:{identity}%] [evalue:{evalue:.2e}]\n")
        f.write(str(bowhead_seq.seq) + "\n")
    print(f"saved! ({identity}%, source: {hit_id})")

### Pairwise alignment for both types of species pairs (cetaceans vs non-cetaceans) ###

print("\n" + "=" * 60)
print("Computing pairwise identity for each species pair")
print("=" * 60)

results = []

for gene in ALL_GENES:
    # figure out if this is a cancer or control gene
    if gene in CANCER_GENES:
        gene_type = "cancer"
    else:
        gene_type = "control"

    ###  Cetacean pair (blue whale vs bowhead whale) 
    blue_file = os.path.join(RAW_FOLDER, f"{gene}_blue_whale.fasta")
    bowhead_file = os.path.join(RAW_FOLDER, f"{gene}_bowhead.fasta")
    blue_seq = read_sequence(blue_file)
    bowhead_seq = read_sequence(bowhead_file)

    # if both sequences exist, run pairwise_identity() and print out percent similarity
    if blue_seq and bowhead_seq:
        cetacean_id = pairwise_identity(blue_seq, bowhead_seq)
        print(f"  {gene:<10} Cetacean: {cetacean_id}%")
    else:
        cetacean_id = "N/A"
        print(f"  {gene:<10} Cetacean: N/A (missing sequence)")

    ### Suidae pair: wild boar vs common warthog
    boar_file = os.path.join(RAW_FOLDER, f"{gene}_wild_boar.fasta")
    warthog_file = os.path.join(RAW_FOLDER, f"{gene}_warthog.fasta")
    boar_seq = read_sequence(boar_file)
    warthog_seq = read_sequence(warthog_file)

    # if both sequences exist, run pairwise_identity() and print out percent similarity
    if boar_seq and warthog_seq:
        suidae_id = pairwise_identity(boar_seq, warthog_seq)
        print(f"  {'':10} Suidae:     {suidae_id}%")
    else:
        suidae_id = "N/A"
        print(f"  {'':10} Suidae:     N/A (missing sequence)")

    # store this gene's results and information
    results.append({
        "Gene": gene,
        "Gene_Type": gene_type,
        "Cetacean_Identity": cetacean_id,
        "Suid_Identity": suidae_id,
    })

# save results to a CSV file
csv_file = os.path.join(RESULTS_FOLDER, "comparison.csv")
with open(csv_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Gene", "Gene_Type", "Cetacean_Identity", "Suid_Identity"])
    writer.writeheader()
    for row in results:
        writer.writerow(row)

print(f"\nResults saved to {csv_file}")
print("=" * 65)
print("NEXT STEP: Run statistical_comparison.py")
print("=" * 65)
