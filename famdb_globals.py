import logging
import re

LOGGER = logging.getLogger(__name__)

GROUP_FAMILIES = "Families"
GROUP_LOOKUP_BYNAME = "Lookup/ByName"
GROUP_LOOKUP_BYACC = "Lookup/ByAccession"
GROUP_LOOKUP_BYSTAGE = "Lookup/ByStage"
GROUP_NODES = "Taxonomy/Nodes"
GROUP_TAXANAMES = "Partitions"

# DF####### or DF########## or DR####### or DR##########
dfam_acc_pat = re.compile("^(D[FR])([0-9]{2})([0-9]{2})([0-9]{2})[0-9]{3,6}$")

# The current version of the file format
FILE_VERSION = "1.0"

# The version of the famdb python package
GENERATOR_VERSION = "1.0.2"

LEAF_LINK = "leaf_link:"
ROOT_LINK = "root_link:"

# TODO: Finalize the location of the dfam export utilities
REPBASE_FILE = "Libraries/RMRB_spec_to_tax.json"

# TODO: make command-line options to customize these
DESCRIPTION = (
    "Dfam - A database of transposable element (TE) sequence alignments and HMMs."
)

COPYRIGHT_TEXT = """Dfam - A database of transposable element (TE) sequence alignments and HMMs
Copyright (C) %s The Dfam consortium.

Release: Dfam_%s
Date   : %s

This database is free; you can redistribute it and/or modify it
as you wish, under the terms of the CC0 1.0 license, a
'no copyright' license:

The Dfam consortium has dedicated the work to the public domain, waiving
all rights to the work worldwide under copyright law, including all related
and neighboring rights, to the extent allowed by law.

You can copy, modify, distribute and perform the work, even for
commercial purposes, all without asking permission.
See Other Information below.

Other Information

o In no way are the patent or trademark rights of any person affected by
  CC0, nor are the rights that other persons may have in the work or in how
  the work is used, such as publicity or privacy rights.
o Makes no warranties about the work, and disclaims liability for all uses of the
  work, to the fullest extent permitted by applicable law.
o When using or citing the work, you should not imply endorsement by the Dfam consortium.

You may also obtain a copy of the CC0 license here:
http://creativecommons.org/publicdomain/zero/1.0/legalcode
"""

FILE_DESCRIPTION = f"""This is famdb.py version {GENERATOR_VERSION}.

example commands, including the most commonly used options:

  famdb.py [-i DB_DIR] info
    Prints information about the file including database name and date.

  famdb.py [-i DB_DIR] names 'mus' | head
    Prints taxonomy nodes that include 'mus', and the corresponding IDs.
    The IDs and names are stored in the FamDB file, and are based
    on the NCBI taxonomy database (https://www.ncbi.nlm.nih.gov/taxonomy).

  famdb.py [-i DB_DIR] lineage -ad 'Homo sapiens'
  famdb.py [-i DB_DIR] lineage -ad --format totals 9606
    Prints a taxonomic tree including the given clade and optionally ancestors
    and/or descendants, with the number of repeats indicated at each level of
    the hierarchy. With the 'totals' format, prints the number of matching
    ancestral and lineage-specific entries.

  famdb.py [-i DB_DIR] family --format fasta_acc MIR3
    Exports a single family from the database in one of several formats.

  famdb.py [-i DB_DIR] families -f embl_meta -ad --curated 'Drosophila melanogaster'
  famdb.py [-i DB_DIR] families -f hmm -ad --curated --class LTR 7227
    Searches and exports multiple families from the database, in one of several formats.

"""

FAMILY_FORMATS_EPILOG = """
Supported formats:
  * 'summary'     : (default) A human-readable summary format. Currently includes
                    accession, name, classification, and length.

  * 'hmm'         : The family's HMM, including some additional metadata such as
                    species and RepeatMasker classification.
  * 'hmm_species' : Same as 'hmm', but with a species-specific TH line extracted
                    into the GA/TC/NC values. This format is only useful for the
                    families command when querying within a species for which such
                    thresholds have been determined.

  * 'fasta_name'  : FASTA, with the following header format:
                    >MIR @Mammalia [S:40,60,65]
  * 'fasta_acc'   : FASTA, with the following header format:
                    >DF0000001.4 @Mammalia [S:40,60,65]

  * 'embl'        : EMBL, including all metadata and the consensus sequence.
  * 'embl_meta'   : Same as 'embl', but with only the metadata included.
  * 'embl_seq'    : Same as 'embl', but with only the sequences included.
"""

MISSING_FILE = """
\tTaxon in Partition %s, Partition File Not Found
\t      
\tThis partition is available for download from Dfam and may
\tbe installed in the %s directory.  Please see
\t%s for more details."""

HELP_URL = "http://www.dfam.org/home"

# Soundex codes
SOUNDEX_LOOKUP = {
    "A": 0,
    "E": 0,
    "I": 0,
    "O": 0,
    "U": 0,
    "Y": 0,
    "B": 1,
    "F": 1,
    "P": 1,
    "V": 1,
    "C": 2,
    "G": 2,
    "J": 2,
    "K": 2,
    "Q": 2,
    "S": 2,
    "X": 2,
    "Z": 2,
    "D": 3,
    "T": 3,
    "L": 4,
    "M": 5,
    "N": 5,
    "R": 6,
    "H": None,
    "W": None,
}
