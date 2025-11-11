"""
Hackbio Internship
Team: Glycine

# Task: DNA to Protein Conversion
# Author: Opeoluwa Shodipe
# Github: https://github.com/Opeoluwa-Shodipe/Single-Cell-RNA-Sequencing/tree/main/Stage%20one
# Linkedin: https://www.linkedin.com/in/sopeoluwa/
"""

from Bio.Seq import Seq

def translate_dna_to_protein(dna, start_at_atg=False, to_stop=False):
    """
    Translate a DNA sequence into its protein sequence using BioPython.

    Args:
        dna (str): DNA sequence (A, T, C, G).
        start_at_atg (bool): If True, translation starts at the first 'ATG'.
        to_stop (bool): If True, translation stops at the first stop codon.

    Returns:
        str: Amino acid sequence (using standard translation table).
    """
    dna_clean = dna.upper().replace(' ', '').replace('\n', '')
    if start_at_atg:
        idx = dna_clean.find('ATG')
        if idx >= 0:
            dna_clean = dna_clean[idx:]
    seq_obj = Seq(dna_clean)
    protein = seq_obj.translate(to_stop=to_stop)
    return str(protein)

# Standard translation
prot = translate_dna_to_protein("ACCGTAACCATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print(prot)  

# Start at ATG and stop at the first stop codon
prot2 = translate_dna_to_protein("ACCGTAACCATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 
                                 start_at_atg=True, to_stop=True)
print(prot2)  
