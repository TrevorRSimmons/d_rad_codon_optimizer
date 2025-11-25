# optimizer.py

# Codon usage table sourced from http://lowelab.ucsc.edu/GtRNAdb/Dein_radi/Dein_radi-summary-codon.html
codon_usage_table = {
    'A': {'GCT': 0.88, 'GCC': 5.97, 'GCG': 4.74, 'GCA': 0.63},
    'G': {'GGT': 0.73, 'GGC': 5.85, 'GGG': 2.08, 'GGA': 0.54},
    'P': {'CCT': 0.52, 'CCC': 2.74, 'CCG': 2.5, 'CCA': 0.3},
    'T': {'ACT': 0.36, 'ACC': 3.47, 'ACG': 1.73, 'ACA': 0.21},
    'V': {'GTT': 0.36, 'GTC': 2.31, 'GTG': 4.84, 'GTA': 0.23},
    'S': {'TCT': 0.2, 'TCC': 0.76, 'TCG': 1.16, 'TCA': 0.17, 'AGT': 0.43, 'AGC': 2.38},
    'R': {'CGT': 0.68, 'CGC': 4.1, 'CGG': 2.06, 'CGA': 0.23, 'AGG': 0.21, 'AGA': 0.11},
    'L': {'CTT': 0.49, 'CTC': 3.43, 'CTG': 6.75, 'CTA': 0.13, 'TTG': 0.73, 'TTA': 0.07},
    'F': {'TTT': 1.05, 'TTC': 2.07},
    'N': {'AAT': 0.35, 'AAC': 2.03},
    'K': {'AAG': 1.94, 'AAA': 0.78},
    'D': {'GAT': 0.53, 'GAC': 4.52},
    'E': {'GAG': 3.22, 'GAA': 2.52},
    'H': {'CAT': 0.33, 'CAC': 1.76},
    'Q': {'CAG': 3.4, 'CAA': 0.71},
    'I': {'ATT': 0.99, 'ATC': 2.19, 'ATA': 0.08},
    'M': {'ATG': 1.77},
    'Y': {'TAT': 0.3, 'TAC': 1.92},
    'C': {'TGT': 0.1, 'TGC': 0.55},
    'W': {'TGG': 1.38},
    '*': {'TAG': 0.03, 'TAA': 0.08, 'TGA': 0.21},
}

def validate_dna_sequence(sequence):
    valid_nucleotides = set("ATGCatgc")
    return all(nucleotide in valid_nucleotides for nucleotide in sequence)

def validate_aa_sequence(sequence):
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    return all(aa in valid_amino_acids for aa in sequence)

def translate_codon(codon):
    translation_table = {
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GGT': 'G', 'GGC': 'G', 'GGG': 'G', 'GGA': 'G',
        'CCT': 'P', 'CCC': 'P', 'CCG': 'P', 'CCA': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACG': 'T', 'ACA': 'T',
        'GTT': 'V', 'GTC': 'V', 'GTG': 'V', 'GTA': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCG': 'S', 'TCA': 'S', 'AGT': 'S', 'AGC': 'S',
        'CGT': 'R', 'CGC': 'R', 'CGG': 'R', 'CGA': 'R', 'AGG': 'R', 'AGA': 'R',
        'CTT': 'L', 'CTC': 'L', 'CTG': 'L', 'CTA': 'L', 'TTG': 'L', 'TTA': 'L',
        'TTT': 'F', 'TTC': 'F',
        'AAT': 'N', 'AAC': 'N',
        'AAG': 'K', 'AAA': 'K',
        'GAT': 'D', 'GAC': 'D',
        'GAG': 'E', 'GAA': 'E',
        'CAT': 'H', 'CAC': 'H',
        'CAG': 'Q', 'CAA': 'Q',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
        'ATG': 'M',
        'TAT': 'Y', 'TAC': 'Y',
        'TGT': 'C', 'TGC': 'C',
        'TGG': 'W',
        'TAG': '*', 'TAA': '*', 'TGA': '*',
    }
    return translation_table.get(codon, '')

def optimize_codons_dna(dna_sequence, codon_usage_table):
    if not validate_dna_sequence(dna_sequence):
        raise ValueError("Invalid DNA sequence. Only A, T, G, and C allowed.")

    optimized_sequence = ""

    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i:i+3]
        if len(codon) != 3:
            continue

        amino_acid = translate_codon(codon.upper())

        if amino_acid == '*':
            optimized_sequence += codon.upper()
        elif amino_acid in codon_usage_table:
            optimal_codon = max(codon_usage_table[amino_acid],
                                key=codon_usage_table[amino_acid].get)
            optimized_sequence += optimal_codon
        else:
            optimized_sequence += codon.upper()

    return optimized_sequence

def translate_aa_to_dna(aa_sequence):
    dna_sequence = ''
    for aa in aa_sequence:
        if aa in codon_usage_table:
            if aa == '*':
                dna_sequence += 'TAA'
            else:
                optimal_codon = max(codon_usage_table[aa],
                                    key=codon_usage_table[aa].get)
                dna_sequence += optimal_codon
    return dna_sequence

def optimize_codons_aa(aa_sequence, codon_usage_table):
    if not validate_aa_sequence(aa_sequence):
        raise ValueError("Invalid amino acid sequence.")
    dna_sequence = translate_aa_to_dna(aa_sequence)
    return optimize_codons_dna(dna_sequence, codon_usage_table)

def optimize_sequence(seq: str, input_type: str) -> str:
    if input_type == "aa":
        return optimize_codons_aa(seq, codon_usage_table)
    else:
        return optimize_codons_dna(seq, codon_usage_table)
