import subprocess
from Bio.Blast import NCBIXML
from Bio import SeqIO
import os
import shutil

# manually add the BLAST+ bin directory to PATH
blast_bin_dir = r"C:\Program Files\NCBI\blast-2.16.0+\bin"
os.environ['PATH'] += os.pathsep + blast_bin_dir

# Was having issues with the program finding my path variables, so i used shutil to force the executable by: 
blastp_exe = shutil.which('blastp')
if blastp_exe is None:
    # fallback to full path if blastp is still not found
    blastp_exe = r"C:\Program Files\NCBI\blast-2.16.0+\bin\blastp.exe"
    if not os.path.exists(blastp_exe):
        raise FileNotFoundError("blastp executable not found. Please ensure BLAST+ is installed and the path is correctly set.")


human_fasta = './human.fa'
mouse_db = './mouse_db'  # the BLAST database created from mouse.fa using in the same directory as mouse.fa: makeblastdb -in mouse.fa -dbtype prot -out mouse_db
output_file = 'blast_results.txt'


with open(output_file, 'w') as out_handle:
    for human_record in SeqIO.parse(human_fasta, 'fasta'):
        query_temp = 'temp_query.fasta'
        SeqIO.write(human_record, query_temp, 'fasta')

        blastp_cmd = [
            blastp_exe,
            '-query', query_temp,
            '-db', mouse_db,
            '-evalue', '0.001',
            '-outfmt', '5',  # XML format
            '-out', 'temp_result.xml',
            '-matrix', 'BLOSUM62'
        ]

        result = subprocess.run(blastp_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


        if result.returncode != 0:
            print(f"Error running BLASTP for sequence {human_record.id}:\n{result.stderr}")
            out_handle.write(f"Human ID: {human_record.id}\n")
            out_handle.write("Error running BLASTP.\n")
            out_handle.write("\n" + "-"*60 + "\n\n")
            os.remove(query_temp)
            continue


        with open('temp_result.xml') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records, None)


            if blast_record and blast_record.alignments:
                top_alignment = blast_record.alignments[0]
                top_hsp = top_alignment.hsps[0]

                human_id = human_record.id
                mouse_id = top_alignment.hit_def.split()[0]  # mouse sequence ID
                e_value = top_hsp.expect
                bitscore = top_hsp.bits
                alignment_query = top_hsp.query
                alignment_match = top_hsp.match
                alignment_sbjct = top_hsp.sbjct

                out_handle.write(f"Human ID: {human_id}\n")
                out_handle.write(f"Mouse ID: {mouse_id}\n")
                out_handle.write(f"E-value: {e_value}\n")
                out_handle.write(f"Bitscore: {bitscore}\n")
                out_handle.write("Alignment:\n")
                out_handle.write(f"Query:  {alignment_query}\n")
                out_handle.write(f"Match:  {alignment_match}\n")
                out_handle.write(f"Sbjct:  {alignment_sbjct}\n")
                out_handle.write("\n" + "-"*60 + "\n\n")
            else:
                out_handle.write(f"Human ID: {human_record.id}\n")
                out_handle.write("No hits found.\n")
                out_handle.write("\n" + "-"*60 + "\n\n")

        # clean up the temp files.
        os.remove(query_temp)
        os.remove('temp_result.xml')


# Explanation of Choices:

# (i) Choice of BLAST Program: 'blastp'
# We are comparing protein sequences from 'human.fa' against a protein database 'mouse_db' created from 'mouse.fa', so we use 'blastp'.
# It identifies homologous proteins based on amino acid sequence similarity.

# (ii) Choice of Substitution Matrix: 'BLOSUM62'
# 'BLOSUM62' is a substitution matrix suitable for detecting medium-range sequence similarities.
# It provides a good balance between sensitivity and specificity, making it ideal for comparing sequences between species like human and mouse.

# (iii) Choice of Parameters:
# The 'blastp_exe' variable is used to specify the path to the BLASTP executable.
# - '-query temp_query.fasta': Specifies the input query file containing the human protein sequences.
#   Each sequence from 'human.fa' is written to 'temp_query.fasta' and used as the query for the BLAST search.
# - '-db mouse_db': Specifies the BLAST database to search against.
#   'mouse_db' is the protein database created from 'mouse.fa', containing mouse protein sequences.
# - '-evalue 0.001': Sets a stringent E-value threshold to report only highly significant matches.
#   A lower E-value (0.001) reduces the chance of false positives, ensuring that only statistically significant alignments are considered.
# - '-outfmt 5': Specifies the output format as XML.
#   XML format (option 5) is structured and easily parsed using Biopython's NCBIXML parser.
# - '-matrix BLOSUM62': Specifies the substitution matrix to use during alignment.
#   Reinforces the choice of 'BLOSUM62' for scoring.
# The default parameters for other options (e.g., gap penalties, word size) are used, which are generally suitable for most protein BLAST searches.