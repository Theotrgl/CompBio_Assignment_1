codon_dict = {
    'UUU': 'Phe (F)', 'UUC': 'Phe (F)', 'UUA': 'Leu (L)', 'UUG': 'Leu (L)',
    'CUU': 'Leu (L)', 'CUC': 'Leu (L)', 'CUA': 'Leu (L)', 'CUG': 'Leu (L)',
    'AUU': 'Ile (I)', 'AUC': 'Ile (I)', 'AUA': 'Ile (I)', 'AUG': 'Met (M)',
    'GUU': 'Val (V)', 'GUC': 'Val (V)', 'GUA': 'Val (V)', 'GUG': 'Val (V)',
    'UCU': 'Ser (S)', 'UCC': 'Ser (S)', 'UCA': 'Ser (S)', 'UCG': 'Ser (S)',
    'CCU': 'Pro (P)', 'CCC': 'Pro (P)', 'CCA': 'Pro (P)', 'CCG': 'Pro (P)',
    'ACU': 'Thr (T)', 'ACC': 'Thr (T)', 'ACA': 'Thr (T)', 'ACG': 'Thr (T)',
    'GCU': 'Ala (A)', 'GCC': 'Ala (A)', 'GCA': 'Ala (A)', 'GCG': 'Ala (A)',
    'UAU': 'Tyr (Y)', 'UAC': 'Tyr (Y)', 'UAA': 'Stop (X)', 'UAG': 'Stop (X)',
    'CAU': 'His (H)', 'CAC': 'His (H)', 'CAA': 'Gln (Q)', 'CAG': 'Gln (Q)',
    'AAU': 'Asn (N)', 'AAC': 'Asn (N)', 'AAA': 'Lys (K)', 'AAG': 'Lys (K)',
    'GAU': 'Asp (D)', 'GAC': 'Asp (D)', 'GAA': 'Glu (E)', 'GAG': 'Glu (E)',
    'UGU': 'Cys (C)', 'UGC': 'Cys (C)', 'UGA': 'Stop (X)', 'UGG': 'Trp (W)',
    'CGU': 'Arg (R)', 'CGC': 'Arg (R)', 'CGA': 'Arg (R)', 'CGG': 'Arg (R)',
    'AGU': 'Ser (S)', 'AGC': 'Ser (S)', 'AGA': 'Arg (R)', 'AGG': 'Arg (R)',
    'GGU': 'Gly (G)', 'GGC': 'Gly (G)', 'GGA': 'Gly (G)', 'GGG': 'Gly (G)'
}

dna_to_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
dna_to_mrna = {'A': 'U', 'T': 'U', 'C': 'G', 'G': 'C'}

def no_1(dna_sequence):
    seqm = dna_sequence.upper()
    complement_sequence = ''.join([dna_to_complement[base] for base in seqm])
    mrna_sequence = ''.join([dna_to_mrna[base] for base in complement_sequence])

    amino_acid_sequence = []
    i = 0
    while i < len(mrna_sequence) - 2:
        codon = mrna_sequence[i:i+3]
        amino_acid = codon_dict.get(codon, 'Unknown')
        amino_acid_sequence.append(amino_acid)
        i += 3
    
    print("Input DNA =", seqm)
    print("Complement =", complement_sequence)
    print("mRNA =", "".join([dna_to_mrna[base] for base in complement_sequence]))
    print("Aminoacid =", " - ".join(amino_acid_sequence))
    
    
no_1("acttaagca")

#%% 
reverse_codon_dict = {
    'N': ['AAU', 'AAC'],
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'C': ['UGU', 'UGC'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['UUU', 'UUC'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'K': ['AAA', 'AAG'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'M': ['AUG'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    'X': ['UAA', 'UAG', 'UGA'],
}

def no_2(amino_acid_chain):
    rna_sequence = ""

    for amino_acid in amino_acid_chain:
        if amino_acid in reverse_codon_dict:
            rna_sequence += reverse_codon_dict[amino_acid][0]

    print("Input Amino-Acid =", amino_acid_chain)
    print("mRNA =", rna_sequence)

    frequency_dict = {}
    i = 0
    while i < len(rna_sequence) - 2:
        codon = rna_sequence[i:i+3]
        frequency_dict[codon] = frequency_dict.get(codon, 0) + 1
        i += 3

    for codon, frequency in frequency_dict.items():
        print(f"{codon} = {frequency}")

no_2("N-A-N")






