"""
blast_compare.py --------------------------------------------------------------------------------------

Author: Paul Antezana 

This program serves to compare the blue whale and bowhead whale gene sequences.

This either uses pairwise alignment (faster, no internet needed) if the sequences are downloaded. 
This works using Biopython which can directly align them and report percent identity

If the sequences are not downloaded, this program utlizes NCBI BLAST using an API which directly
posts to the NCBI's BLAST server through the requests library. This gives us timeout control, something
Biopython's NCBI.WWW lacks in providing. By BLASTing the sequence against NT, we then scan the top 100
hits for any bowhead whale matches.

BEFORE RUNNING:##########
  run fetch_sequences.py
#########################
"""

# import necessary Python modules
import os
import re
import time
import requests
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIXML
from io import StringIO

# build file paths relative to this script's location so it works from any directory
script_dir = os.path.dirname(os.path.abspath(__file__))
raw_folder = os.path.join(script_dir, "../data/raw") + "/"
results_folder = os.path.join(script_dir, "../results") + "/"

# create the results folder if it doesn't already exist yet
os.makedirs(results_folder, exist_ok=True)

# a list of the current genes we're interested in
genes = ["CIRBP", "RPA2", "PCNA", "ERCC1", "TP53"]

# url for NCBI's online BLAST tool
BLAST_URL = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"

# this variable determines how long the program will wait for NCBI to finish a given BLAST job before giving up
BLAST_TIMEOUT = 300

def read_sequence(filepath):
    """
    Purpose: Read a FASTA file and return its DNA sequence as a plain string.

    Parameters:
        filepath (str): the path to the FASTA file to read

    Returns:
        str: the DNA sequence if the file exists
        None: if the file was not found
    """
    if not os.path.exists(filepath):
        return None
    record = SeqIO.read(filepath, "fasta")
    return str(record.seq)

def pairwise_identity(seq1, seq2):
    """
    Purpose: Align two DNA sequences and calculate how similar they are as a percentage.

    Parameters:
        seq1 (str): the first DNA sequence (e.g. bowhead)
        seq2 (str): the second DNA sequence (e.g. blue whale)

    Returns:
        float: percent identity between the two sequences, rounded to 2 decimal places
    """
    aligner = PairwiseAligner()
    # local alignment handles partial overlaps pretty much like how BLAST does
    aligner.mode = "local"  
    best = next(iter(aligner.align(seq1, seq2)))
    matches = sum(a == b for a, b in zip(*best))
    return round((matches / best.length) * 100, 2)

def submit_blast_job(sequence):
    """
    Purpose: Submit a DNA sequence to NCBI's online BLAST server and get back a job ID.

    Parameters:
        sequence (str): the DNA sequence to search with (blue whale sequence)

    Returns:
        str: the Request ID (RID) for the submitted job if successful
        None: if the submission failed
    """
    # send the sequence to NCBI's BLAST server and get back a job ID (RID)
    params = {
        "CMD":          "Put",
        "PROGRAM":      "blastn",
        "DATABASE":     "nt",
        "QUERY":        sequence,

        # request to access the top 100 hits. this increases the chance of finding bowhead
        "HITLIST_SIZE": "100",   
        "FORMAT_TYPE":  "XML",
        "EXPECT":       "0.001",
    }
    try:
        response = requests.post(BLAST_URL, data=params, timeout=30)
    except requests.exceptions.RequestException:
        return None

    # NCBI embeds the job ID in the response to something like: RID = ABCDE11111
    match = re.search(r"RID = (\w+)", response.text)
    if match:
        return match.group(1)
    return None

def poll_blast_results(request_id):
    """
    Purpose: Repeatedly check NCBI database for the results of a submitted BLAST job until they are ready.

    Parameters:
        request_id (str): the Request ID returned by submit_blast_job()

    Returns:
        str: the full BLAST results as an XML string if successful
        None: if results timed out or failed
    """
    # keep checking NCBI every 20 seconds until the results are ready or we time out
    start = time.time()

    while time.time() - start < BLAST_TIMEOUT:
        # NCBI asks you to wait at least ~15 seconds between checks
        time.sleep(20)  
        elapsed = int(time.time() - start)
        print(f"    checking ({elapsed}s)...", end=" ", flush=True)

        try:
            response = requests.get(
                BLAST_URL,
                params={"CMD": "Get", "RID": request_id, "FORMAT_TYPE": "XML"},
                timeout=60
            )
        except requests.exceptions.RequestException as e:
            print(f"network error: {e}")
            continue

        text = response.text

        if "Status=FAILED" in text or "Status=UNKNOWN" in text:
            print("NCBI reported failure.")
            return None

        # when results are ready, NCBI returns XML starting with <BlastOutput
        if "<BlastOutput" in text:
            print("done!")
            return text

    # if BLAST takes too long (i.e. surpasses the BLAST_TIMEOUT variable), we print out this message and return None.
    print("timed out.")
    return None


def find_bowhead_identity(xml_text):
    """
    Purpose: Parse BLAST XML results and find the best matching hit from the bowhead whale.

    Parameters:
        xml_text (str): the raw XML string returned by poll_blast_results()

    Returns:
        float: percent identity of the best bowhead hit, or None if no bowhead hit was found
    """
    blast_record = NCBIXML.read(StringIO(xml_text))

    for alignment in blast_record.alignments:
        title = alignment.title.lower()
        # this code runs if we found something belonging to the bowhead whale
        if "balaena mysticetus" in title or "bowhead" in title:
            hsp = alignment.hsps[0]
            identity = (hsp.identities / hsp.align_length) * 100
            return round(identity, 2)
        
    # no bowhead hit found in the top 100 results
    return None  

def blast_vs_bowhead(sequence):
    """
    Purpose: Run the full BLAST pipeline to find the bowhead whale match for a given sequence.

    Parameters:
        sequence (str): the blue whale DNA sequence to use as the query

    Returns:
        float: percent identity to the best bowhead hit
        None: no match was found
    """
    # run the full BLAST pipeline:
    # --> submit a job --> wait for results --> find bowhead hit from results
    rid = submit_blast_job(sequence)
    if rid is None:
        return None

    xml_text = poll_blast_results(rid)
    if xml_text is None:
        return None

    return find_bowhead_identity(xml_text)

print("### Step 2: Comparing gene sequences (bowhead vs blue whale) ###\n")

# empty dictionary for results
all_results = {}

# this runs for all the genes in our genes list, one at a time
for gene in genes:
    blue_file    = f"{raw_folder}{gene}_blue_whale.fasta"
    bowhead_file = f"{raw_folder}{gene}_bowhead.fasta"
    blue_seq    = read_sequence(blue_file)
    bowhead_seq = read_sequence(bowhead_file)
    print(f"{gene}...", end=" ", flush=True)

    # if the blue and bowhead whale sequences exist, then they are both on disk. if this is true, 
    # use pairwise_alignment (it's faster and doesn't require internet)
    if blue_seq and bowhead_seq:
        identity = pairwise_identity(bowhead_seq, blue_seq)
        print(f"{identity}%  (pairwise alignment)")
        all_results[gene] = identity

    # if ONLY the blue whale sequence exists, use NCBI BLAST to find the bowhead match
    elif blue_seq:
        identity = blast_vs_bowhead(blue_seq)
        if identity is not None:
            print(f"{identity}%  (NCBI BLAST)")
        else:
            print("No bowhead hit found in top 100 results.")
        all_results[gene] = identity

    else:
        # neither sequence file exists, therefore, skip this gene
        print("No sequence file found --- skipping...")
        all_results[gene] = None

# create a .csv file to store our result data in
output_file = f"{results_folder}blast_results.csv"

# write down results into the .csv file
with open(output_file, "w") as f:
    f.write("Gene,Percent_Identity\n")
    for gene, identity in all_results.items():
        value = str(identity) if identity is not None else "N/A"
        f.write(f"{gene},{value}\n")

# if the program ran without major errors, then these print statements should work. 
# this serves as a way to verify the program worked.
print(f"\nAll done! Results saved to {output_file}")
print("#########################################################")
print("NEXT STEPS: Run analyze_results.py to generate the chart.")
print("#########################################################")
