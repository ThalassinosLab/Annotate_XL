# Monoisotopic element masses
MASS_DICT = {
    "H": 1.007825032,
    "C": 12.0096,
    "O": 15.99491462,
    "N": 14.00307401,
    "S": 31.97207069
}

# Monoisotopic amino acid masses as they appear in peptide backbone
AMINO_MONO_MASS = {
    'G': 57.021464,
    'A': 71.037114,
    'S': 87.032028,
    'P': 97.052764,
    'V': 99.068414,
    'M': 131.040485,
    'X': 147.0354,  # single oxidised methionine sulfoxide.
    'T': 101.04768,
    'H': 137.058912,
    'C': 160.030644,  # with carbamidomethylation (57.02146+103.009184)
    'F': 147.068414,
    'L': 113.084064,
    'R': 156.101111,
    'I': 113.084064,
    'N': 114.042927,
    'D': 115.026943,
    'Q': 128.058578,
    'K': 128.094963,
    'E': 129.042593,
    'W': 186.079313,
    'Y': 163.063329
}

# Immonium ion masses
IMMON_MASS = {
    'G': 30.03438,
    'A': 44.05003,
    'S': 60.04494,
    'P': 70.06568,
    'V': 72.08133,
    'M': 104.0534,
    'X': 120.0483,  # single oxidised methionine sulfoxide.
    'T': 74.06059,
    'H': 110.0718,
    'C': 133.0436,  # with carbamidomethylation (57.02146+103.009184)
    'F': 120.0813,
    'L': 86.09698,
    'R': 129.114,
    'I': 86.09698,
    'N': 87.05584,
    'D': 88.03986,
    'Q': 101.0715,
    'K': 84.08136,
    'E': 102.0555,
    'W': 159.0922,
    'Y': 136.0762
}

# Unique diagnostic ions formed by BS3/DSS crosslinker
DIAG_IONS = {
    'DI_1': 139.08,
    'DI_2': 222.15
}