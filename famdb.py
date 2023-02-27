#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    famdb.py, version 0.4.3
    Usage: famdb.py [-h] [-l LOG_LEVEL] [-i FILE] command ...

    Queries or modifies the contents of a famdb file. For more detailed help
    and information about program options, run `famdb.py --help` or
    `famdb.py <command> --help`.

    This program can also be used as a module. It provides classes and methods
    for working with FamDB files, which contain Transposable Element (TE)
    families and associated taxonomy data.

    # Classes
        Family: Metadata and model of a TE family.
        FamDB: HDF5-based format for storing Family objects.

SEE ALSO:
    Dfam: http://www.dfam.org

AUTHOR(S):
    Jeb Rosen <jeb.rosen@systemsbiology.org>

LICENSE:
    This code may be used in accordance with the Creative Commons
    Zero ("CC0") public domain dedication:
    https://creativecommons.org/publicdomain/zero/1.0/

DISCLAIMER:
    This software is provided ``AS IS'' and any express or implied
    warranties, including, but not limited to, the implied warranties of
    merchantability and fitness for a particular purpose, are disclaimed.
    In no event shall the authors or the Dfam consortium members be
    liable for any direct, indirect, incidental, special, exemplary, or
    consequential damages (including, but not limited to, procurement of
    substitute goods or services; loss of use, data, or profits; or
    business interruption) however caused and on any theory of liability,
    whether in contract, strict liability, or tort (including negligence
    or otherwise) arising in any way out of the use of this software, even
    if advised of the possibility of such damage.
"""

import argparse
import collections
import datetime
import json
import logging
import os
import re
import sys
import textwrap
import time

import h5py
import numpy


LOGGER = logging.getLogger(__name__)


# Soundex codes
SOUNDEX_LOOKUP = {
    'A': 0, 'E': 0, 'I': 0, 'O': 0, 'U': 0, 'Y': 0,
    'B': 1, 'F': 1, 'P': 1, 'V': 1,
    'C': 2, 'G': 2, 'J': 2, 'K': 2, 'Q': 2, 'S': 2, 'X': 2, 'Z': 2,
    'D': 3, 'T': 3,
    'L': 4,
    'M': 5, 'N': 5,
    'R': 6,
    'H': None, 'W': None,
}

def soundex(word):
    """
    Converts 'word' according to American Soundex[1].

    This is used for "sounds like" types of searches.

    [1]: https://en.wikipedia.org/wiki/Soundex#American_Soundex
    """

    codes = [SOUNDEX_LOOKUP[ch] for ch in word.upper() if ch in SOUNDEX_LOOKUP]

    # Start at the second code
    i = 1

    # Drop identical sounds and H and W
    while i < len(codes):
        code = codes[i]
        prev = codes[i-1]

        if code is None:
            # Drop H and W
            del codes[i]
        elif code == prev:
            # Drop adjacent identical sounds
            del codes[i]
        else:
            i += 1

    # Keep the first letter
    coding = word[0]

    # Keep codes, except for the first or vowels
    codes_rest = filter(lambda c: c > 0, codes[1:])

    # Append stringified remaining numbers
    for code in codes_rest:
        coding += str(code)

    # Pad to 3 digits
    while len(coding) < 4:
        coding += '0'

    # Truncate to 3 digits
    return coding[:4]

def sounds_like(first, second):
    """
    Returns true if the string 'first' "sounds like" 'second'.

    The comparison is currently implemented by running both strings through the
    soundex algorithm and checking if the soundex values are equal.
    """
    soundex_first = soundex(first)
    soundex_second = soundex(second)

    return soundex_first == soundex_second

def sanitize_name(name):
    """
    Returns the "sanitized" version of the given 'name'.
    This must be kept in sync with Dfam's algorithm.
    """
    name = re.sub(r"[\s\,\_]+", "_", name)
    name = re.sub(r"[\(\)\<\>\']+", "", name)
    return name

class Family:  # pylint: disable=too-many-instance-attributes
    """A Transposable Element family, made up of metadata and a model."""

    FamilyField = collections.namedtuple("FamilyField", ["name", "type"])

    # Known metadata fields
    META_FIELDS = [
        # Core required family metadata
        FamilyField("name", str),
        FamilyField("accession", str),
        FamilyField("version", int),
        FamilyField("consensus", str),
        FamilyField("length", int),

        # Optional family metadata
        FamilyField("title", str),
        FamilyField("author", str),
        FamilyField("description", str),
        FamilyField("classification", str),
        FamilyField("classification_note", str),
        FamilyField("search_stages", str),
        FamilyField("buffer_stages", str),
        FamilyField("clades", list),
        FamilyField("date_created", str),
        FamilyField("date_modified", str),
        FamilyField("repeat_type", str),
        FamilyField("repeat_subtype", str),
        FamilyField("features", str),
        FamilyField("coding_sequences", str),
        FamilyField("aliases", str),
        FamilyField("citations", str),
        FamilyField("refineable", bool),
        FamilyField("target_site_cons", str),
        # TODO: add source_assembly, source_method?

        # Metadata available when a model is present
        FamilyField("model", str),
        FamilyField("max_length", int),
        FamilyField("is_model_masked", bool),
        FamilyField("seed_count", int),
        FamilyField("build_method", str),
        FamilyField("search_method", str),
        FamilyField("taxa_thresholds", str),
        FamilyField("general_cutoff", float),
    ]

    # Metadata lookup by field name
    META_LOOKUP = {field.name: field for field in META_FIELDS}

    @staticmethod
    def type_for(name):
        """Returns the expected data type for the attribute 'name'."""
        return Family.META_LOOKUP[name].type

    def __getattr__(self, name):
        if name not in Family.META_LOOKUP:
            raise AttributeError("Unknown Family metadata attribute '{}'".format(name))

    # Data is converted on setting, so that consumers can rely on the correct types
    def __setattr__(self, name, value):
        if name not in Family.META_LOOKUP:
            raise AttributeError("Unknown Family metadata attribute '{}'".format(name))

        expected_type = self.type_for(name)
        if value is not None and not isinstance(value, expected_type):
            try:
                value = expected_type(value)
            except Exception as exc:
                raise TypeError("Incompatible type for '{}'. Expected '{}', got '{}'".format(
                    name, expected_type, type(value))) from exc
        super().__setattr__(name, value)

    def accession_with_optional_version(self):
        """
        Returns the accession of 'self', with '.version' appended if the version is known.
        """

        acc = self.accession
        if self.version is not None:
            acc += "." + str(self.version)
        return acc

    # A useful string representation for debugging, but not much else
    def __str__(self):
        return "%s.%s '%s': %s len=%d" % (self.accession, self.version,
                                          self.name, self.classification, self.length or -1)

    def to_dfam_hmm(self, famdb, species=None, include_class_in_name=False, require_general_threshold=False):  # pylint: disable=too-many-locals,too-many-branches
        """
        Converts 'self' to Dfam-style HMM format.
        'famdb' is used for lookups in the taxonomy database (id -> name).

        If 'species' (a taxonomy id) is given, the assembly-specific GA/TC/NC
        thresholds will be used instead of the threshold that was in the HMM
        (usually a generic or strictest threshold).
        """
        if self.model is None:
            return None

        out = ""

        # Appends to 'out':
        # "TAG   Text"
        #
        # Or if wrap=True and 'text' has multiple lines:
        # "TAG   Line 1"
        # "TAG   Line 2"
        def append(tag, text, wrap=False):
            nonlocal out
            if not text:
                return

            prefix = "%-6s" % tag
            text = str(text)
            if wrap:
                text = textwrap.fill(text, width=72)
            out += textwrap.indent(text, prefix)
            out += "\n"

        # TODO: Compare to e.g. finditer(). This does a lot of unnecessary
        # allocation since most of model_lines are appended verbatim.
        model_lines = self.model.split("\n")

        i = 0
        for i, line in enumerate(model_lines):
            if line.startswith("HMMER3"):
                out += line + "\n"

                name = self.name or self.accession
                if include_class_in_name:
                    rm_class = self.repeat_type
                    if self.repeat_subtype:
                        rm_class += "/" + self.repeat_subtype
                    name = name + "#" + rm_class

                append("NAME", name)
                append("ACC", self.accession_with_optional_version())
                append("DESC", self.title)
            elif any(map(line.startswith, ["NAME", "ACC", "DESC"])):
                # Correct version of this line was output already
                pass
            elif line.startswith("CKSUM"):
                out += line + "\n"
                break
            else:
                out += line + "\n"

        th_lines = []
        species_hmm_ga = None
        species_hmm_tc = None
        species_hmm_nc = None
        if self.taxa_thresholds:
            for threshold in self.taxa_thresholds.split("\n"):
                parts = threshold.split(",")
                tax_id = int(parts[0])
                (hmm_ga, hmm_tc, hmm_nc, hmm_fdr) = map(float, parts[1:])

                tax_name = famdb.get_taxon_name(tax_id, 'scientific name')
                if tax_id == species:
                    species_hmm_ga, species_hmm_tc, species_hmm_nc = hmm_ga, hmm_tc, hmm_nc
                th_lines += ["TaxId:%d; TaxName:%s; GA:%.2f; TC:%.2f; NC:%.2f; fdr:%.3f;" % (
                    tax_id, tax_name, hmm_ga, hmm_tc, hmm_nc, hmm_fdr)]

        if species is None:
            if self.general_cutoff:
                species_hmm_ga = species_hmm_tc = species_hmm_nc = self.general_cutoff

        if species_hmm_ga:
            append("GA", "%.2f;" % species_hmm_ga)
            append("TC", "%.2f;" % species_hmm_tc)
            append("NC", "%.2f;" % species_hmm_nc)
        elif require_general_threshold:
            LOGGER.debug("missing general threshold for " + self.accession)
            return None

        for th_line in th_lines:
            append("TH", th_line)

        if self.build_method:
            append("BM", self.build_method)
        if self.search_method:
            append("SM", self.search_method)

        append("CT", (self.classification and self.classification.replace("root;", "")))

        for clade_id in self.clades:
            tax_name = famdb.get_sanitized_name(clade_id)
            append("MS", "TaxId:%d TaxName:%s" % (clade_id, tax_name))

        append("CC", self.description, True)
        append("CC", "RepeatMasker Annotations:")
        append("CC", "     Type: %s" % (self.repeat_type or ""))
        append("CC", "     SubType: %s" % (self.repeat_subtype or ""))

        species_names = [famdb.get_sanitized_name(c) for c in self.clades]
        append("CC", "     Species: %s" % ", ".join(species_names))

        append("CC", "     SearchStages: %s" % (self.search_stages or ""))
        append("CC", "     BufferStages: %s" % (self.buffer_stages or ""))

        if self.refineable:
            append("CC", "     Refineable")

        # Append all remaining lines unchanged
        out += "\n".join(model_lines[i+1:])

        return out

    __COMPLEMENT_TABLE = str.maketrans(
        "ACGTRYWSKMNXBDHV",
        "TGCAYRSWMKNXVHDB"
    )

    def to_fasta(
            self,
            famdb,
            use_accession=False,
            include_class_in_name=False,
            do_reverse_complement=False,
            buffer=None
    ):
        """Converts 'self' to FASTA format."""
        sequence = self.consensus
        if sequence is None:
            return None
        sequence = sequence.upper()

        if use_accession:
            identifier = self.accession_with_optional_version()
        else:
            identifier = self.name or self.accession

        if buffer:
            if buffer is True:
                # range-less specification: leave identifier unchanged, and use
                # the whole sequence as the buffer
                buffer = [1, len(sequence)]
            else:
                # range specification: append _START_END to the identifier
                identifier += "_%d_%d" % (buffer[0], buffer[1])

            sequence = sequence[buffer[0]-1:buffer[1]]
            identifier = identifier + "#buffer"

        if do_reverse_complement:
            sequence = sequence.translate(self.__COMPLEMENT_TABLE)
            sequence = sequence[::-1]

        if include_class_in_name and not buffer:
            rm_class = self.repeat_type
            if self.repeat_subtype:
                rm_class += "/" + self.repeat_subtype
            identifier = identifier + "#" + rm_class

        header = ">" + identifier

        if do_reverse_complement:
            header += " (anti)"

        if use_accession and self.name:
            header += " name=" + self.name

        for clade_id in self.clades:
            clade_name = famdb.get_sanitized_name(clade_id)
            header += " @" + clade_name

        if self.search_stages:
            header += " [S:%s]" % self.search_stages

        out = header + "\n"

        i = 0
        while i < len(sequence):
            out += sequence[i:i+60] + "\n"
            i += 60

        return out

    def to_embl(self, famdb, include_meta=True, include_seq=True):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
        """Converts 'self' to EMBL format."""

        if include_seq and self.consensus is None:
            # Skip families without consensus sequences, if sequences were required.
            # metadata-only formats will still include families without a consensus sequence.
            return None

        sequence = self.consensus or ""

        out = ""

        # Appends to 'out':
        # "TAG  Text"
        #
        # Or if wrap=True and 'text' has multiple lines:
        # "TAG  Line 1"
        # "TAG  Line 2"
        def append(tag, text, wrap=False):
            nonlocal out
            if not text:
                return

            prefix = "%-5s" % tag
            if wrap:
                text = textwrap.fill(str(text), width=72)
            out += textwrap.indent(str(text), prefix)
            out += "\n"

        # Appends to 'out':
        # "FT                   line 1"
        # "FT                   line 2"
        def append_featuredata(text):
            nonlocal out
            prefix = "FT                   "
            if text:
                out += textwrap.indent(textwrap.fill(str(text), width=72), prefix)
                out += "\n"

        id_line = self.accession
        if self.version is not None:
            id_line += "; SV " + str(self.version)

        append("ID", "%s; linear; DNA; STD; UNC; %d BP." % (id_line, len(sequence)))
        append("NM", self.name)
        out += "XX\n"
        append("AC", self.accession + ';')
        out += "XX\n"
        append("DE", self.title, True)
        out += "XX\n"

        if include_meta:
            if self.aliases:
                for alias_line in self.aliases.splitlines():
                    [db_id, db_link] = map(str.strip, alias_line.split(":"))
                    if db_id == "Repbase":
                        append("DR", "Repbase; %s." % db_link)
                        out += "XX\n"

            if self.repeat_type == "LTR":
                append("KW", "Long terminal repeat of retrovirus-like element; %s." % self.name)
            else:
                append("KW", "%s/%s." % (self.repeat_type or "", self.repeat_subtype or ""))
            out += "XX\n"

            for clade_id in self.clades:
                lineage = famdb.get_lineage_path(clade_id)
                if lineage[0] == "root":
                    lineage = lineage[1:]

                if len(lineage) > 0:
                    append("OS", lineage[-1])
                    append("OC", "; ".join(lineage[:-1]) + ".", True)
            out += "XX\n"

            if self.citations:
                citations = json.loads(self.citations)
                citations.sort(key=lambda c: c["order_added"])
                for cit in citations:
                    append("RN", "[%d] (bases 1 to %d)" % (cit["order_added"], self.length))
                    append("RA", cit["authors"], True)
                    append("RT", cit["title"], True)
                    append("RL", cit["journal"])
                    out += "XX\n"

            append("CC", self.description, True)
            out += "CC\n"
            append("CC", "RepeatMasker Annotations:")
            append("CC", "     Type: %s" % (self.repeat_type or ""))
            append("CC", "     SubType: %s" % (self.repeat_subtype or ""))

            species_names = [famdb.get_sanitized_name(c)
                             for c in self.clades]
            append("CC", "     Species: %s" % ", ".join(species_names))

            append("CC", "     SearchStages: %s" % (self.search_stages or ""))
            append("CC", "     BufferStages: %s" % (self.buffer_stages or ""))
            if self.refineable:
                append("CC", "     Refineable")

            if self.coding_sequences:
                out += "XX\n"
                append("FH", "Key             Location/Qualifiers")
                out += "FH\n"
                for cds in json.loads(self.coding_sequences):
                    # TODO: sanitize values which might already contain a " in them?

                    append("FT", "CDS             %d..%d" % (cds["cds_start"], cds["cds_end"]))
                    append_featuredata('/product="%s"' % cds["product"])
                    append_featuredata('/number=%s' % cds["exon_count"])
                    append_featuredata('/note="%s"' % cds["description"])
                    append_featuredata('/translation="%s"' % cds["translation"])


            out += "XX\n"

        if include_seq:
            sequence = sequence.lower()
            i = 0
            counts = {"a": 0, "c": 0, "g": 0, "t": 0, "other": 0}
            for char in sequence:
                if char not in counts:
                    char = "other"
                counts[char] += 1

            append("SQ", "Sequence %d BP; %d A; %d C; %d G; %d T; %d other;" % (
                len(sequence), counts["a"], counts["c"], counts["g"], counts["t"],
                counts["other"]))

            while i < len(sequence):
                chunk = sequence[i:i+60]
                i += 60

                j = 0
                line = ""
                while j < len(chunk):
                    line += chunk[j:j + 10] + " "
                    j += 10

                out += "     %-66s %d\n" % (line, min(i, len(sequence)))

        out += "//\n"

        return out

    @staticmethod
    def read_embl_families(filename, lookup, header_cb=None):
        """
        Iterates over Family objects from the .embl file 'filename'. The format
        should match the output format of to_embl(), but this is not thoroughly
        tested.

        'lookup' should be a dictionary of Species names (in the EMBL file) to
        taxonomy IDs.

        If specified, 'header_cb' will be invoked with the contents of the
        header text at the top of the file before the iteration is complete.

        TODO: This mechanism is a bit awkward and should perhaps be reworked.
        """

        def set_family_code(family, code, value):
            """
            Sets an attribute on 'family' based on the EMBL line starting with 'code'.
            For codes corresponding to list attributes, values are appended.
            """
            if code == "ID":
                match = re.match(r'(\S*)', value)
                acc = match.group(1)
                acc = acc.rstrip(";")
                family.accession = acc
            elif code == "NM":
                family.name = value
            elif code == "DE":
                family.description = value
            elif code == "CC":
                # TODO: Consider only recognizing these after seeing "RepeatMasker Annotations"

                matches = re.match(r'\s*Type:\s*(\S+)', value)
                if matches:
                    family.repeat_type = matches.group(1).strip()

                matches = re.match(r'\s*SubType:\s*(\S+)', value)
                if matches:
                    family.repeat_subtype = matches.group(1).strip()

                matches = re.search(r'Species:\s*(.+)', value)
                if matches:
                    for spec in matches.group(1).split(","):
                        name = spec.strip().lower()
                        if name:
                            tax_id = lookup.get(name)
                            if tax_id:
                                family.clades += [tax_id]
                            else:
                                LOGGER.warning("Could not find taxon for '%s'", name)

                matches = re.search(r'SearchStages:\s*(\S+)', value)
                if matches:
                    family.search_stages = matches.group(1).strip()

                matches = re.search(r'BufferStages:\s*(\S+)', value)
                if matches:
                    family.buffer_stages = matches.group(1).strip()

                matches = re.search(r'Refineable', value)
                if matches:
                    family.refineable = True

        header = ""
        family = None
        in_header = True
        in_metadata = False

        with open(filename) as file:
            for line in file:
                if family is None:
                    # ID indicates start of metadata
                    if line.startswith("ID"):
                        family = Family()
                        family.clades = []
                        in_header = False
                        in_metadata = True
                    elif in_header:
                        matches = re.match(r"(CC)?\s*(.*)", line)
                        if line.startswith("XX"):
                            in_header = False
                        elif matches:
                            header_line = matches.group(2).rstrip('*').strip()
                            header += header_line + "\n"
                        else:
                            header += line

                if family is not None:
                    if in_metadata:
                        # SQ line indicates start of sequence
                        if line.startswith("SQ"):
                            in_metadata = False
                            family.consensus = ""

                        # Continuing metadata
                        else:
                            split = line.rstrip("\n").split(None, maxsplit=1)
                            if len(split) > 1:
                                code = split[0].strip()
                                value = split[1].strip()
                                set_family_code(family, code, value)

                    # '//' line indicates end of the sequence area
                    elif line.startswith("//"):
                        family.length = len(family.consensus)

                        yield family
                        family = None

                    # Part of the sequence area
                    else:
                        family.consensus += re.sub(r'[^A-Za-z]', '', line)

        if header_cb:
            header_cb(header)


# The current version of the file format
FILE_VERSION = "0.5"

class FamDB:
    """Transposable Element Family and taxonomy database."""

    dtype_str = h5py.special_dtype(vlen=str)

    GROUP_FAMILIES = "Families"
    GROUP_FAMILIES_BYNAME = "Families/ByName"
    GROUP_FAMILIES_BYACC = "Families/ByAccession"
    GROUP_FAMILIES_BYSTAGE = "Families/ByStage"
    GROUP_NODES = "Taxonomy/Nodes"

    # DF####### or DF########## or DR####### or DR##########
    dfam_acc_pat = re.compile("^(D[FR])([0-9]{2})([0-9]{2})[0-9]{3,6}$")

    def __init__(self, filename, mode="r"):
        if mode == "r":
            reading = True

            # If we definitely will not be writing to the file, optimistically assume
            # nobody else is writing to it and disable file locking. File locking can
            # be a bit flaky, especially on NFS, and is unnecessary unless there is
            # a parallel writer (which is unlikely for famdb files).
            os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

        elif mode == "r+":
            reading = True
        elif mode == "w":
            reading = False
        else:
            raise ValueError("Invalid file mode. Expected 'r' or 'r+' or 'w', got '{}'".format(mode))


        self.filename = filename
        self.file = h5py.File(filename, mode)
        self.mode = mode

        try:
            if reading and self.file.attrs["version"] != FILE_VERSION:
                raise Exception("File version is {}, but this is version {}".format(
                    self.file.attrs["version"], FILE_VERSION,
                ))
        except:
            # This 'except' catches both "version" missing from attrs, or the
            # value not matching if it is present.
            raise Exception("This file cannot be read by this version of famdb.py.")

        self.__lineage_cache = {}

        if self.mode == "w":
            self.seen = {}
            self.added = {'consensus': 0, 'hmm': 0}
            self.__write_metadata()
        elif self.mode == "r+":
            self.seen = {}
            self.seen["name"] = set(self.file[FamDB.GROUP_FAMILIES_BYNAME].keys())
            self.seen["accession"] = set(self.__families_iterator(self.file[FamDB.GROUP_FAMILIES],"Families"))
            self.added = self.get_counts()

        if reading:
            self.names_dump = json.loads(self.file["TaxaNames"][0])

    def __write_metadata(self):
        self.file.attrs["generator"] = "famdb.py v0.4.3"
        self.file.attrs["version"] = FILE_VERSION
        self.file.attrs["created"] = str(datetime.datetime.now())

    def __write_counts(self):
        self.file.attrs["count_consensus"] = self.added['consensus']
        self.file.attrs["count_hmm"] = self.added['hmm']

    def set_db_info(self, name, version, date, desc, copyright_text):
        """Sets database metadata for the current file"""
        self.file.attrs["db_name"] = name
        self.file.attrs["db_version"] = version
        self.file.attrs["db_date"] = date
        self.file.attrs["db_description"] = desc
        self.file.attrs["db_copyright"] = copyright_text

    def get_db_info(self):
        """
        Gets database database metadata for the current file as a dict with keys
        'name', 'version', 'date', 'description', 'copyright'
        """
        if "db_name" not in self.file.attrs:
            return None

        return {
            "name": self.file.attrs["db_name"],
            "version": self.file.attrs["db_version"],
            "date": self.file.attrs["db_date"],
            "description": self.file.attrs["db_description"],
            "copyright": self.file.attrs["db_copyright"],
        }

    def get_metadata(self):
        """
        Gets file metadata for the current file as a dict with keys
        'generator', 'version', 'created'
        """
        return {
            "generator": self.file.attrs["generator"],
            "version": self.file.attrs["version"],
            "created": self.file.attrs["created"],
        }

    def get_counts(self):
        """
        Gets counts of entries in the current file as a dict
        with 'consensus', 'hmm'
        """
        return {
            "consensus": self.file.attrs["count_consensus"],
            "hmm": self.file.attrs["count_hmm"],
        }

    def close(self):
        """Closes this FamDB instance, making further use invalid."""
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __check_unique(self, family, key):
        """Verifies that 'family' is uniquely identified by its value of 'key'."""

        seen = self.seen
        value = getattr(family, key)
        if key not in seen:
            seen[key] = set()

        if value in seen[key]:
            raise Exception("Family is not unique! Already seen {}: {}".format(key, value))

        seen[key].add(value)

    @staticmethod
    def __accession_bin(acc):
        """Maps an accession (Dfam or otherwise) into apropriate bins (groups) in HDF5"""
        dfam_match = FamDB.dfam_acc_pat.match(acc)
        if dfam_match:
            path = FamDB.GROUP_FAMILIES + "/" + dfam_match.group(1) + "/" + \
                   dfam_match.group(2) + "/" + dfam_match.group(3)
        else:
            path = FamDB.GROUP_FAMILIES + "/Aux/" + acc[0:2].lower()
        return path

    @staticmethod
    def __families_iterator(g, prefix=''):
        for key, item in g.items():
            path = '{}/{}'.format(prefix, key)
            if isinstance(item, h5py.Dataset): # test for dataset
                yield (key)
            elif isinstance(item, h5py.Group): # test for group (go down)
                yield from FamDB.__families_iterator(item, path)

    def add_family(self, family):
        """Adds the family described by 'family' to the database."""
        # Verify uniqueness of name and accession.
        # This is important because of the links created to them later.
        if family.name:
            self.__check_unique(family, "name")
        self.__check_unique(family, "accession")

        # Increment counts
        if family.consensus:
            self.added['consensus'] += 1
        if family.model:
            self.added['hmm'] += 1

        # Create the family data
        # In v0.5 we bin the datasets into subgroups to improve performance
        group_path = self.__accession_bin(family.accession)
        dset = self.file.require_group(group_path).create_dataset(family.accession, (0,))

        # Set the family attributes
        for k in Family.META_LOOKUP:
            value = getattr(family, k)
            if value:
                if k == "model":
                    #self.file['hmm_dataset'][self.added['hmm']] = value
                    dset.attrs[k] = value
                else:
                    dset.attrs[k] = value

        # Create links
        if family.name:
            self.file.require_group(FamDB.GROUP_FAMILIES_BYNAME)[family.name] = h5py.SoftLink(group_path + "/" + family.accession)
        # In FamDB format version 0.5 we removed the /Families/ByAccession group as it's redundant
        # (all the data is in Families/<datasets> *and* HDF5 suffers from poor performance when
        # the number of entries in a group exceeds 200-500k.

        for clade_id in family.clades:
            taxon_group = self.file.require_group(FamDB.GROUP_NODES).require_group(str(clade_id))
            families_group = taxon_group.require_group("Families")
            families_group[family.accession] = h5py.SoftLink(group_path + "/" + family.accession)

        def add_stage_link(stage, accession):
            stage_group = self.file.require_group(FamDB.GROUP_FAMILIES_BYSTAGE).require_group(stage.strip())
            if accession not in stage_group:
                stage_group[accession] = h5py.SoftLink(group_path + "/" + accession)

        if family.search_stages:
            for stage in family.search_stages.split(","):
                add_stage_link(stage, family.accession)

        if family.buffer_stages:
            for stage in family.buffer_stages.split(","):
                stage = stage.split("[")[0]
                add_stage_link(stage, family.accession)

        LOGGER.debug("Added family %s (%s)", family.name, family.accession)

    def write_taxonomy(self, tax_db):
        """Writes taxonomy nodes in 'tax_db' to the database."""
        LOGGER.info("Writing taxonomy nodes to database")
        start = time.perf_counter()

        self.names_dump = {}

        count = 0
        for taxon in tax_db.values():
            if taxon.used:
                count += 1

                self.names_dump[taxon.tax_id] = taxon.names

        def store_tree_links(taxon, parent_id):
            group = self.file.require_group(FamDB.GROUP_NODES).require_group(str(taxon.tax_id))
            if parent_id:
                group.create_dataset("Parent", data=[parent_id])

            child_ids = []
            for child in taxon.children:
                if child.used:
                    child_ids += [child.tax_id]
                    store_tree_links(child, taxon.tax_id)

            group.create_dataset("Children", data=child_ids)

        names_data = numpy.array([json.dumps(self.names_dump)])
        names_dset = self.file.create_dataset("TaxaNames", shape=names_data.shape,
                                              dtype=FamDB.dtype_str)
        names_dset[:] = names_data

        LOGGER.info("Writing taxonomy tree")
        # 1 is the "root" taxon
        store_tree_links(tax_db[1], None)

        delta = time.perf_counter() - start
        LOGGER.info("Wrote %d taxonomy nodes in %f", count, delta)

    def finalize(self):
        """Writes some collected metadata, such as counts, to the database"""

        self.__write_counts()

    def has_taxon(self, tax_id):
        """Returns True if 'self' has a taxonomy entry for 'tax_id'"""
        return str(tax_id) in self.file[FamDB.GROUP_NODES]

    def search_taxon_names(self, text, kind=None, search_similar=False):
        """
        Searches 'self' for taxons with a name containing 'text', returning an
        iterator that yields a tuple of (id, is_exact) for each matching node.
        Each id is returned at most once, and if any of its names are an exact
        match the whole node is treated as an exact match.

        If 'similar' is True, names that sound similar will also be considered
        eligible.

        A list of strings may be passed as 'kind' to restrict what kinds of
        names will be searched.
        """

        text = text.lower()

        for tax_id, names in self.names_dump.items():
            matches = False
            exact = False
            for name_cls, name_txt in names:
                name_txt = name_txt.lower()
                if kind is None or kind == name_cls:
                    if text == name_txt:
                        matches = True
                        exact = True
                    elif name_txt.startswith(text + " <"):
                        matches = True
                        exact = True
                    elif text == sanitize_name(name_txt):
                        matches = True
                        exact = True
                    elif text in name_txt:
                        matches = True
                    elif search_similar and sounds_like(text, name_txt):
                        matches = True

            if matches:
                yield [int(tax_id), exact]

    def resolve_species(self, term, kind=None, search_similar=False):
        """
        Resolves 'term' as a species or clade in 'self'. If 'term' is a number,
        it is a taxon id. Otherwise, it will be searched for in 'self' in the
        name fields of all taxa. A list of strings may be passed as 'kind' to
        restrict what kinds of names will be searched.

        If 'search_similar' is True, a "sounds like" search will be tried
        first. If it is False, a "sounds like" search will still be performed

        if no results were found.

        This function returns a list of tuples (taxon_id, is_exact) that match
        the query. The list will be empty if no matches were found.
        """

        # Try as a number
        try:
            tax_id = int(term)
            if self.has_taxon(tax_id):
                return [[tax_id, True]]

            return []
        except ValueError:
            pass

        # Perform a search by name, splitting between exact and inexact matches for sorting
        exact = []
        inexact = []
        for tax_id, is_exact in self.search_taxon_names(term, kind, search_similar):
            if is_exact:
                exact += [tax_id]
            else:
                inexact += [tax_id]

        # Combine back into one list, with exact matches first
        results = [[tax_id, True] for tax_id in exact]
        for tax_id in inexact:
            results += [[tax_id, False]]

        if len(results) == 0 and not search_similar:
            # Try a sounds-like search (currently soundex)
            similar_results = self.resolve_species(term, kind, True)
            if similar_results:
                print("No results were found for that name, but some names sound similar:",
                      file=sys.stderr)
                for tax_id, _ in similar_results:
                    names = self.get_taxon_names(tax_id)
                    print(tax_id, ", ".join(["{1}".format(*n) for n in names]), file=sys.stderr)

        return results

    def resolve_one_species(self, term, kind=None):
        """
        Resolves 'term' in 'dbfile' as a taxon id or search term unambiguously.
        Parameters are as in the 'resolve_species' method.
        Returns None if not exactly one result is found,
        and prints details to the screen.
        """

        results = self.resolve_species(term, kind)

        # Check for a single exact match first, to any field
        exact_matches = []
        for nid, is_exact in results:
            if is_exact:
                exact_matches += [nid]
        if len(exact_matches) == 1:
            return exact_matches[0]

        if len(results) == 0:
            print("No species found for search term '{}'".format(term), file=sys.stderr)
        elif len(results) == 1:
            return results[0][0]
        else:
            print("""Ambiguous search term '{}' (found {} results, {} exact).
Please use a more specific name or taxa ID, which can be looked
up with the 'names' command."""
                  .format(term, len(results), len(exact_matches)), file=sys.stderr)

        return None

    def get_taxon_names(self, tax_id):
        """
        Returns a list of [name_class, name_value] of the taxon given by 'tax_id'.
        """

        return self.names_dump[str(tax_id)]

    def get_taxon_name(self, tax_id, kind='scientific name'):
        """
        Returns the first name of the given 'kind' for the taxon given by 'tax_id',
        or None if no such name was found.
        """

        names = self.names_dump[str(tax_id)]
        for name in names:
            if name[0] == kind:
                return name[1]

        return None

    def get_sanitized_name(self, tax_id):
        """
        Returns the "sanitized name" of tax_id, which is the sanitized version
        of the scientific name.
        """

        name = self.get_taxon_name(tax_id, 'scientific name')
        if name:
            name = sanitize_name(name)
        return name

    def get_families_for_taxon(self, tax_id):
        """Returns a list of the accessions for each family directly associated with 'tax_id'."""
        group = self.file[FamDB.GROUP_NODES][str(tax_id)].get("Families")
        if group:
            return group.keys()
        else:
            return []

    def get_lineage(self, tax_id, **kwargs):
        """
        Returns the lineage of 'tax_id'. Recognized kwargs: 'descendants' to include
        descendant taxa, 'ancestors' to include ancestor taxa.
        IDs are returned as a nested list, for example
        [ 1, [ 2, [3, [4]], [5], [6, [7]] ] ]
        where '2' may have been the passed-in 'tax_id'.
        """

        group_nodes = self.file[FamDB.GROUP_NODES]

        if kwargs.get("descendants"):
            def descendants_of(tax_id):
                descendants = [int(tax_id)]
                for child in group_nodes[str(tax_id)]["Children"]:
                    descendants += [descendants_of(child)]
                return descendants
            tree = descendants_of(tax_id)
        else:
            tree = [tax_id]

        if kwargs.get("ancestors"):
            while tax_id:
                node = group_nodes[str(tax_id)]
                if "Parent" in node:
                    tax_id = node["Parent"][0]
                    tree = [tax_id, tree]
                else:
                    tax_id = None

        return tree

    def get_lineage_path(self, tax_id, cache=True):
        """
        Returns a list of strings encoding the lineage for 'tax_id'.
        """

        if cache and tax_id in self.__lineage_cache:
            return self.__lineage_cache[tax_id]

        tree = self.get_lineage(tax_id, ancestors=True)
        lineage = []

        while tree:
            node = tree[0]
            tree = tree[1] if len(tree) > 1 else None

            tax_name = self.get_taxon_name(node, 'scientific name')
            lineage += [tax_name]

        if cache:
            self.__lineage_cache[tax_id] = lineage

        return lineage

    @staticmethod
    def __filter_name(family, name):
        """Returns True if the family's name begins with 'name'."""

        if family.attrs.get("name"):
            if family.attrs["name"].lower().startswith(name):
                return True

        return False

    def __filter_stages(self, accession, stages):
        """Returns True if the family belongs to a search or buffer stage in 'stages'."""
        for stage in stages:
            grp = self.file[FamDB.GROUP_FAMILIES_BYSTAGE].get(stage)
            if grp and accession in grp:
                return True

        return False

    @staticmethod
    def __filter_search_stages(family, stages):
        """Returns True if the family belongs to a search stage in 'stages'."""
        if family.attrs.get("search_stages"):
            sstages = (ss.strip() for ss in family.attrs["search_stages"].split(","))
            for family_ss in sstages:
                if family_ss in stages:
                    return True

        return False

    @staticmethod
    def __filter_repeat_type(family, rtype):
        """
        Returns True if the family's RepeatMasker Type plus SubType
        (e.g. "DNA/CMC-EnSpm") starts with 'rtype'.
        """
        if family.attrs.get("repeat_type"):
            full_type = family.attrs["repeat_type"]
            if family.attrs.get("repeat_subtype"):
                full_type = full_type + "/" + family.attrs["repeat_subtype"]

            if full_type.lower().startswith(rtype):
                return True

        return False

    @staticmethod
    def __filter_curated(accession, curated):
        """
        Returns True if the family's curatedness is the same as 'curated'. In
        other words, 'curated=True' includes only curated familes and
        'curated=False' includes only uncurated families.

        Families are currently assumed to be curated unless their name is of the
        form DR<7 digits>.

        TODO: perhaps this should be a dedicated 'curated' boolean field on Family
        """

        is_curated = not (
            accession.startswith("DR") and
            len(accession) == 9 and
            all((c >= "0" and c <= "9" for c in accession[2:])))

        return is_curated == curated

    def get_accessions_filtered(self, **kwargs):
        """
        Returns an iterator that yields accessions for the given search terms.

        Filters are specified in kwargs:
            tax_id: int
            ancestors: boolean, default False
            descendants: boolean, default False
                If none of (tax_id, ancestors, descendants) are
                specified, *all* families will be checked.
            curated_only = boolean
            uncurated_only = boolean
            stage = int
            is_hmm = boolean
            repeat_type = string (prefix)
            name = string (prefix)
                If any of stage, repeat_type, or name are
                omitted (or None), they will not be used to filter.

                The behavior of 'stage' depends on 'is_hmm': if is_hmm is True,
                stage must match in SearchStages (a match in BufferStages is not
                enough).
        """

        if not ("tax_id" in kwargs or "ancestors" in kwargs or "descendants" in kwargs):
            tax_id = 1
            ancestors = True
            descendants = True
        else:
            tax_id = kwargs["tax_id"]
            ancestors = kwargs.get("ancestors") or False
            descendants = kwargs.get("descendants") or False

        # Define family filters (logically ANDed together)
        filters = []

        if kwargs.get("curated_only"):
            filters += [lambda a, f: self.__filter_curated(a, True)]
        if kwargs.get("uncurated_only"):
            filters += [lambda a, f: self.__filter_curated(a, False)]

        filter_stage = kwargs.get("stage")
        filter_stages = None
        if filter_stage:
            if filter_stage == 80:
                # "stage 80" = "all stages", so skip filtering
                pass
            elif filter_stage == 95:
                # "stage 95" = this specific stage list:
                filter_stages = ["35", "50", "55", "60", "65", "70", "75"]
                filters += [lambda a, f: self.__filter_stages(a, filter_stages)]
            else:
                filter_stages = [str(filter_stage)]
                filters += [lambda a, f: self.__filter_stages(a, filter_stages)]

        # HMM only: add a search stage filter to "un-list" families that were
        # allowed through only because they match in buffer stage
        if kwargs.get("is_hmm") and filter_stages:
            filters += [lambda a, f: self.__filter_search_stages(f(), filter_stages)]

        filter_repeat_type = kwargs.get("repeat_type")
        if filter_repeat_type:
            filter_repeat_type = filter_repeat_type.lower()
            filters += [lambda a, f: self.__filter_repeat_type(f(), filter_repeat_type)]

        filter_name = kwargs.get("name")
        if filter_name:
            filter_name = filter_name.lower()
            filters += [lambda a, f: self.__filter_name(f(), filter_name)]

        # Recursive iterator flattener
        def walk_tree(tree):
            """Returns all elements in 'tree' with all levels flattened."""
            if hasattr(tree, "__iter__"):
                for elem in tree:
                    yield from walk_tree(elem)
            else:
                yield tree

        seen = set()

        def iterate_accs():
            # special case: Searching the whole database in a specific
            # stage only is a common usage pattern in RepeatMasker.
            # When searching the whole database instead of a species,
            # the number of accessions to read through is shorter
            # when going off of only the stage indexes.
            if tax_id == 1 and descendants and filter_stages \
                    and not filter_repeat_type and not filter_name:
                for stage in filter_stages:
                    grp = self.file[FamDB.GROUP_FAMILIES_BYSTAGE].get(stage)
                    if grp:
                        yield from grp.keys()

            # special case: Searching the whole database, going directly via
            # Families/ is faster than repeatedly traversing the tree
            elif tax_id == 1 and descendants:
                #yield from self.file[FamDB.GROUP_FAMILIES_BYACC].keys()
                for name in self.__families_iterator(self.file[FamDB.GROUP_FAMILIES],"Families"):
                    yield name
            else:
                lineage = self.get_lineage(tax_id, ancestors=ancestors, descendants=descendants)
                for node in walk_tree(lineage):
                    yield from self.get_families_for_taxon(node)


        for accession in iterate_accs():
            if accession in seen:
                continue
            seen.add(accession)

            cached_family = None
            def family_getter():
                nonlocal cached_family
                if not cached_family:
                    path = self.__accession_bin(accession)
                    cached_family = self.file[path].get(accession)
                return cached_family

            match = True
            for filt in filters:
                if not filt(accession, family_getter):
                    match = False
            if match:
                yield accession

    def get_family_names(self):
        """Returns a list of names of families in the database."""
        return sorted(self.file[FamDB.GROUP_FAMILIES_BYNAME].keys(), key=str.lower)

    # Currently I do not believe this is used by anyone?
    def get_family_accessions(self):
        """Returns a list of accessions for families in the database."""
        return sorted(self.__families_iterator(self.file[FamDB.GROUP_FAMILIES],"Families"), key=str.lower)

    @staticmethod
    def __get_family(entry):
        if not entry:
            return None

        family = Family()

        # Read the family attributes and data
        for k in entry.attrs:
            value = entry.attrs[k]
            setattr(family, k, value)

        return family

    def get_family_by_accession(self, accession):
        """Returns the family with the given accession."""
        path = self.__accession_bin(accession)
        if path in self.file:
            entry = self.file[path].get(accession)
            return self.__get_family(entry)
        return None

    def get_family_by_name(self, name):
        """Returns the family with the given name."""
        # TODO: This will also suffer the performance issues seen with
        #       other groups that exceed 200-500k entries in a single group
        #       at some point.  This needs to be refactored to scale appropriately.
        entry = self.file["Families/ByName"].get(name)
        return self.__get_family(entry)


# Command-line utilities

def command_info(args):
    """The 'info' command displays some of the stored metadata."""

    db_info = args.file.get_db_info()
    counts = args.file.get_counts()
    f_info = args.file.get_metadata()

    print("""\
File: {}
FamDB Generator: {}
FamDB Format Version: {}
FamDB Creation Date: {}

Database: {}
Version: {}
Date: {}

{}

Total consensus sequences: {}
Total HMMs: {}
""".format(
    os.path.realpath(args.file.filename),
    f_info["generator"], f_info["version"],
    f_info["created"],
    db_info["name"], db_info["version"],
    db_info["date"], db_info["description"],
    counts["consensus"], counts["hmm"]))

def command_names(args):
    """The 'names' command displays all names of all taxa that match the search term."""

    entries = []
    for tax_id, is_exact in args.file.resolve_species(args.term):
        names = args.file.get_taxon_names(tax_id)
        entries += [[tax_id, is_exact, names]]

    if args.format == "pretty":
        prev_exact = None
        for (tax_id, is_exact, names) in entries:
            if is_exact != prev_exact:
                if is_exact:
                    print("Exact Matches\n=============")
                else:
                    if prev_exact:
                        print()
                    print("Non-exact Matches\n=================")
                prev_exact = is_exact

            print(tax_id, ", ".join(["{1} ({0})".format(*n) for n in names]))
    elif args.format == "json":
        obj = []
        for (tax_id, _, names) in entries:
            names_obj = [{"kind": name[0], "value": name[1]} for name in names]
            obj += [{"id": tax_id, "names": names_obj}]
        print(json.dumps(obj))
    else:
        raise ValueError("Unimplemented names format: %s" % args.format)


def print_lineage_tree(file, tree, gutter_self, gutter_children):
    """Pretty-prints a lineage tree with box drawing characters."""
    if not tree:
        return

    tax_id = tree[0]
    children = tree[1:]
    name = file.get_taxon_name(tax_id, 'scientific name')
    count = len(file.get_families_for_taxon(tax_id))
    print("{}{} {} [{}]".format(gutter_self, tax_id, name, count))

    # All but the last child need a downward-pointing line that will link up
    # to the next child, so this is split into two cases
    if len(children) > 1:
        for child in children[:-1]:
            print_lineage_tree(file, child, gutter_children + "├─", gutter_children + "│ ")

    if children:
        print_lineage_tree(file, children[-1], gutter_children + "└─", gutter_children + "  ")


def print_lineage_semicolons(file, tree, parent_name, starting_at):
    """
    Prints a lineage tree as a flat list of semicolon-delimited names.

    In order to print the correct lineage string, the available tree must
    be "complete" even if ancestors were not specified to build up the
    string starting from "root". 'starting_at' specifies the first taxa
    (in the descending direction) to actually be output.
    """
    if not tree:
        return

    tax_id = tree[0]
    children = tree[1:]
    name = file.get_taxon_name(tax_id, 'scientific name')
    if parent_name:
        name = parent_name + ";" + name

    if starting_at == tax_id:
        starting_at = None

    if not starting_at:
        count = len(file.get_families_for_taxon(tax_id))
        print("{}: {} [{}]".format(tax_id, name, count))

    for child in children:
        print_lineage_semicolons(file, child, name, starting_at)

def get_lineage_totals(file, tree, target_id, seen=None):
    """
    Recursively calculates the total number of families
    on ancestors and descendants of 'target_id' in the given 'tree'.

    'seen' is required to track families that are present on multiple
    lineages due to horizontal transfer and ensure each family
    is only counted one time, either as an ancestor or a descendant.
    """
    if not seen:
        seen = set()

    tax_id = tree[0]
    children = tree[1:]
    accessions = file.get_families_for_taxon(tax_id)

    count_here = 0
    for acc in accessions:
        if acc not in seen:
            seen.add(acc)
            count_here += 1

    if target_id == tax_id:
        target_id = None

    counts = [0, 0]
    for child in children:
        anc, desc = get_lineage_totals(file, child, target_id, seen)
        counts[0] += anc
        counts[1] += desc

    if target_id is None:
        counts[1] += count_here
    else:
        counts[0] += count_here

    return counts

def command_lineage(args):
    """The 'lineage' command outputs ancestors and/or descendants of the given taxon."""

    # TODO: like 'families', filter curated or uncurated (and other filters?)

    target_id = args.file.resolve_one_species(args.term)
    if not target_id:
        return
    tree = args.file.get_lineage(target_id,
                                 descendants=args.descendants,
                                 ancestors=args.ancestors or args.format == "semicolon")

    # TODO: prune branches with 0 total
    if args.format == "pretty":
        print_lineage_tree(args.file, tree, "", "")
    elif args.format == "semicolon":
        print_lineage_semicolons(args.file, tree, "", target_id)
    elif args.format == "totals":
        totals = get_lineage_totals(args.file, tree, target_id)
        print("{} entries in ancestors; {} lineage-specific entries"
              .format(totals[0], totals[1]))
    else:
        raise ValueError("Unimplemented lineage format: %s" % args.format)

def print_families(args, families, header, species=None):
    """
    Prints each family in 'families', optionally with a copyright header. The
    format is determined by 'args.format' and additional data (such as
    taxonomy) is taken from 'args.file'.

    If 'species' is provided and the format is "hmm_species", it is the id of
    the taxa whose species-specific thresholds should be substituted into the
    GA, NC, and TC lines of the HMM.
    """

    # These args are only available with the "families" command. When
    # print_families is called by the "family" command, accessing e.g.
    # args.stage directly raises an AttributeError
    # TODO: consider reworking argument passing to avoid this workaround
    add_reverse_complement = getattr(args, "add_reverse_complement", False)
    include_class_in_name = getattr(args, "include_class_in_name", False)
    require_general_threshold = getattr(args, "require_general_threshold", False)
    stage = getattr(args, "stage", None)

    if header:
        db_info = args.file.get_db_info()
        if db_info:
            copyright_text = db_info["copyright"]
            # Add appropriate comment character to the copyright header lines
            if "hmm" in args.format:
                copyright_text = re.sub("(?m)^", "#   ", copyright_text)
            elif "fasta" in args.format:
                copyright_text = None
            elif "embl" in args.format:
                copyright_text = re.sub("(?m)^", "CC   ", copyright_text)
            if copyright_text:
                print(copyright_text)

    for family in families:
        if args.format == "summary":
            entry = str(family) + "\n"
        elif args.format == "hmm":
            entry = family.to_dfam_hmm(args.file, include_class_in_name=include_class_in_name, require_general_threshold=require_general_threshold)
        elif args.format == "hmm_species":
            entry = family.to_dfam_hmm(args.file, species, include_class_in_name=include_class_in_name, require_general_threshold=require_general_threshold)
        elif args.format == "fasta" or args.format == "fasta_name" or args.format == "fasta_acc":
            use_accession = (args.format == "fasta_acc")

            buffers = []
            if stage and family.buffer_stages:
                for spec in family.buffer_stages.split(","):
                    if "[" in spec:
                        matches = re.match(r'(\d+)\[(\d+)-(\d+)\]', spec.strip())
                        if matches:
                            if stage == int(matches.group(1)):
                                buffers += [[int(matches.group(2)), int(matches.group(3))]]
                        else:
                            LOGGER.warning("Ingored invalid buffer specification: '%s'",
                                           spec.strip())
                    else:
                        buffers += [stage == int(spec)]

            if not buffers:
                buffers += [None]

            entry = ""
            for buffer_spec in buffers:
                entry += family.to_fasta(
                    args.file,
                    use_accession=use_accession,
                    include_class_in_name=include_class_in_name,
                    buffer=buffer_spec
                ) or ""

                if add_reverse_complement:
                    entry += family.to_fasta(args.file,
                                             use_accession=use_accession,
                                             include_class_in_name=include_class_in_name,
                                             do_reverse_complement=True,
                                             buffer=buffer_spec) or ""
        elif args.format == "embl":
            entry = family.to_embl(args.file)
        elif args.format == "embl_meta":
            entry = family.to_embl(args.file, include_meta=True, include_seq=False)
        elif args.format == "embl_seq":
            entry = family.to_embl(args.file, include_meta=False, include_seq=True)
        else:
            raise ValueError("Unimplemented family format: %s" % args.format)

        if entry:
            print(entry, end="")


def command_family(args):
    """The 'family' command outputs a single family by name or accession."""
    family = args.file.get_family_by_accession(args.accession)
    if not family:
        family = args.file.get_family_by_name(args.accession)

    if family:
        print_families(args, [family], False)


def command_families(args):
    """The 'families' command outputs all families associated with the given taxon."""
    target_id = args.file.resolve_one_species(args.term)
    if not target_id:
        return

    families = []

    is_hmm = args.format.startswith("hmm")

    # NB: This is speed-inefficient, because get_accessions_filtered needs to
    # read the whole family data even though we read it again right after.
    # However it is *much* more memory-efficient than loading all the family
    # data at once and then sorting by accession.
    accessions = sorted(args.file.get_accessions_filtered(tax_id=target_id,
                                                          descendants=args.descendants,
                                                          ancestors=args.ancestors,
                                                          curated_only=args.curated,
                                                          uncurated_only=args.uncurated,
                                                          is_hmm=is_hmm,
                                                          stage=args.stage,
                                                          repeat_type=args.repeat_type,
                                                          name=args.name))

    families = map(args.file.get_family_by_accession, accessions)

    print_families(args, families, True, target_id)

def command_append(args):
    """
    The 'append' command reads an EMBL file and appends its entries to an
    existing famdb file.
    """

    lookup = {}
    for tax_id, names in args.file.names_dump.items():
        for name in names:
            if name[0] == "scientific name":
                sanitized_name = sanitize_name(name[1]).lower()
                lookup[sanitized_name] = int(tax_id)

    header = None
    def set_header(val):
        nonlocal header
        header = val

    embl_iter = Family.read_embl_families(args.infile, lookup, set_header)

    seen_accs = args.file.seen["accession"]
    seen_names = args.file.seen["name"]

    for entry in embl_iter:
        acc = entry.accession
        # TODO: This is awkward. The EMBL files being appended may only have an
        # "accession", but that accession may match the *name* of a family
        # already in Dfam. The accession may also match a family already in
        # Dfam, but with a "v" added.
        if acc in seen_accs or acc in seen_names \
            or acc + "v" in seen_accs or acc + "v" in seen_names:
            LOGGER.debug("Ignoring duplicate entry %s", entry.accession)
        else:
            args.file.add_family(entry)

    db_info = args.file.get_db_info()

    if args.name:
        db_info["name"] = args.name
    if args.description:
        db_info["description"] += "\n" + args.description

    db_info["copyright"] += "\n\n" + header

    args.file.set_db_info(
                          db_info["name"], db_info["version"], db_info["date"],
                          db_info["description"], db_info["copyright"]
    )

    # Write the updated counts and metadata
    args.file.finalize()

def main():
    """Parses command-line arguments and runs the requested command."""

    logging.basicConfig()

    parser = argparse.ArgumentParser(description="""This is famdb.py version 0.4.3.

example commands, including the most commonly used options:

  famdb.py [-i FILE] info
    Prints information about the file including database name and date.

  famdb.py [-i FILE] names 'mus' | head
    Prints taxonomy nodes that include 'mus', and the corresponding IDs.
    The IDs and names are stored in the FamDB file, and are based
    on the NCBI taxonomy database (https://www.ncbi.nlm.nih.gov/taxonomy).

  famdb.py [-i FILE] lineage -ad 'Homo sapiens'
  famdb.py [-i FILE] lineage -ad --format totals 9606
    Prints a taxonomic tree including the given clade and optionally ancestors
    and/or descendants, with the number of repeats indicated at each level of
    the hierarchy. With the 'totals' format, prints the number of matching
    ancestral and lineage-specific entries.

  famdb.py [-i FILE] family --format fasta_acc MIR3
    Exports a single family from the database in one of several formats.

  famdb.py [-i FILE] families -f embl_meta -ad --curated 'Drosophila melanogaster'
  famdb.py [-i FILE] families -f hmm -ad --curated --class LTR 7227
    Searches and exports multiple families from the database, in one of several formats.

""", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-l", "--log-level", default="INFO")

    parser.add_argument("-i", "--file", help="specifies the file to query")

    subparsers = parser.add_subparsers(description="""Specifies the kind of query to perform.
For more information on all the possible options for a command, add the --help option after it:
famdb.py families --help
""")

    p_info = subparsers.add_parser("info", description="List general information about the file.")
    p_info.set_defaults(func=command_info)

    p_names = subparsers.add_parser("names", description="List the names and taxonomy identifiers of a clade.")
    p_names.add_argument("-f", "--format", default="pretty", choices=["pretty", "json"],
                         metavar="<format>",
                         help="choose output format. The default is 'pretty'. 'json' is more appropriate for scripts.")
    p_names.add_argument("term", nargs="+", help="search term. Can be an NCBI taxonomy identifier or part of a scientific or common name")
    p_names.set_defaults(func=command_names)

    p_lineage = subparsers.add_parser("lineage", description="List the taxonomy tree including counts of families at each clade.")
    p_lineage.add_argument("-a", "--ancestors", action="store_true",
                           help="include all ancestors of the given clade")
    p_lineage.add_argument("-d", "--descendants", action="store_true",
                           help="include all descendants of the given clade")
    p_lineage.add_argument("-f", "--format", default="pretty", choices=["pretty", "semicolon", "totals"],
                           metavar="<format>",
                           help="choose output format. The default is 'pretty'. 'semicolon' is more appropriate for scripts. 'totals' displays the number of ancestral and lineage-specific families found.")
    p_lineage.add_argument("term", nargs="+", help="search term. Can be an NCBI taxonomy identifier or an unambiguous scientific or common name")
    p_lineage.set_defaults(func=command_lineage)

    family_formats = ["summary", "hmm", "hmm_species", "fasta_name", "fasta_acc", "embl", "embl_meta", "embl_seq"]
    family_formats_epilog = """
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

    p_families = subparsers.add_parser("families", description="Retrieve the families associated \
with a given clade, optionally filtered by additional criteria",
                                       epilog=family_formats_epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    p_families.add_argument("-a", "--ancestors", action="store_true",
                            help="include all ancestors of the given clade")
    p_families.add_argument("-d", "--descendants", action="store_true",
                            help="include all descendants of the given clade")
    p_families.add_argument("--stage", type=int,
                            help="include only families that should be searched in the given stage")
    p_families.add_argument("--class", dest="repeat_type", type=str,
                            help="include only families that have the specified repeat Type/SubType")
    p_families.add_argument("--name", type=str,
                            help="include only families whose name begins with this search term")
    p_families.add_argument("--uncurated", action="store_true",
                            help="include only 'uncurated' families (i.e. named DRXXXXXXX)")
    p_families.add_argument("--curated", action="store_true",
                            help="include only 'curated' families (i.e. not named DRXXXXXXX)")
    p_families.add_argument("-f", "--format", default="summary", choices=family_formats,
                            metavar="<format>", help="choose output format.")
    p_families.add_argument("--add-reverse-complement", action="store_true", help="include a reverse-complemented copy of each matching family; only suppported for fasta formats")
    p_families.add_argument("--include-class-in-name", action="store_true", help="include the RepeatMasker type/subtype after the name (e.g. HERV16#LTR/ERVL); only supported for hmm and fasta formats")
    p_families.add_argument("--require-general-threshold", action="store_true", help="skip families missing general thresholds (and log their accessions at the debug log level)")
    p_families.add_argument("term", nargs="+", help="search term. Can be an NCBI taxonomy identifier or an unambiguous scientific or common name")
    p_families.set_defaults(func=command_families)

    p_family = subparsers.add_parser("family", description="Retrieve details of a single family.",
                                     epilog=family_formats_epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    p_family.add_argument("-f", "--format", default="summary", choices=family_formats,
                          metavar="<format>", help="choose output format.")
    p_family.add_argument("accession", help="the accession of the family to be retrieved")
    p_family.set_defaults(func=command_family)

    p_append = subparsers.add_parser("append")
    p_append.add_argument("infile", help="the name of the input file to be appended")
    p_append.add_argument("--name", help="new name for the database (replaces the existing name)")
    p_append.add_argument("--description", help="additional database description (added to the existing description)")
    p_append.set_defaults(func=command_append)

    args = parser.parse_args()
    logging.getLogger().setLevel(getattr(logging, args.log_level.upper()))

    if "term" in args:
        args.term = " ".join(args.term)

    # For RepeatMasker: Try Libraries/RepeatMaskerLib.h5, if no file was specified
    # in the arguments and that file exists.
    if not args.file:
        # sys.path[0], if non-empty, is initially set to the directory of the
        # originally-invoked script.
        if sys.path[0]:
            default_file = os.path.join(sys.path[0], "Libraries/RepeatMaskerLib.h5")
            if os.path.exists(default_file):
                args.file = default_file

    if args.file:
        try:
            if "func" in args and args.func is command_append:
                mode = "r+"
            else:
                mode = "r"

            args.file = FamDB(args.file, mode)
        except:
            args.file = None

            exc_value = sys.exc_info()[1]
            LOGGER.error("Error reading file: %s", exc_value)
            if LOGGER.getEffectiveLevel() <= logging.DEBUG:
                raise
    else:
        LOGGER.error("Please specify a file to operate on with the -i/--file option.")

    if not args.file:
        return

    if "func" in args:
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # This workaround is from
        # https://docs.python.org/3/library/signal.html#note-on-sigpipe

        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE
