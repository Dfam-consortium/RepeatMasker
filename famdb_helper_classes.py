# coding: utf-8
import collections
import textwrap
import json
import re
import sys


from famdb_globals import LOGGER, LEAF_LINK, ROOT_LINK


class Lineage(list):  # TODO replace exits  with real exception
    """A class to mediate lineages across multiple FamDB files. Contains methods to combine lineages at cross-file break points"""

    def __init__(self, lineage, root, partition_num):
        super().__init__(lineage)
        self.ancestors = False
        self.descendants = False
        self.root = root
        self.partition = partition_num
        links = {LEAF_LINK: {}, ROOT_LINK: None}

        # find and store link elements
        lin_str = lineage.__repr__()
        splits = lin_str.split("'")
        if LEAF_LINK in lin_str:
            self.descendants = True
            if not self.root:
                raise Exception("Leaf Links Found In Non-Root Lineage")
            for i in range(len(splits)):
                if LEAF_LINK in splits[i]:
                    links[LEAF_LINK][i] = str(splits[i].split(":")[1])
        elif ROOT_LINK in lin_str:
            self.ancestors = True
            if self.root:
                raise Exception("Root Links found In Root Lineage")
            links[ROOT_LINK] = {str(lineage[1][0]): str(lineage[1])}

        if self.ancestors and self.descendants:
            raise Exception("Lineage Should Not Contain Root Links And Leaf Links")
        self.links = links
        self.splits = splits

    def __add__(self, other):
        # check to avoid adding root+root or leaf+leaf
        if (self.root and other.root) or (not self.root and not other.root):
            raise Exception("Must Combine Root and Non-Root Lineages")

        # assign lineages
        root_lineage = self if self.root else other
        leaf_lineage = self if not self.root else other
        # load links and split lineage string
        leaf_links = root_lineage.links[LEAF_LINK]
        full_lineage = root_lineage.splits.copy()
        # check each link position in root for linked subtree in leaf
        for position in leaf_links:
            node = leaf_links[position]
            subtree = leaf_lineage.links[ROOT_LINK].get(node)
            if subtree:
                # replace link position with subtree if found
                full_lineage[position] = subtree
        # format splits for json reading
        for i in range(len(full_lineage)):
            if "leaf_link:" in full_lineage[i]:
                full_lineage[i] = f'"{full_lineage[i]}"'
        # join splits back to single string, read as list
        linked_lineage = json.loads("".join(full_lineage))
        # return as new Lineage object, with root True and partion 0
        return Lineage(linked_lineage, root_lineage.root, root_lineage.partition)

    def __iadd__(self, other):
        self.lineage = self + other
        return self.lineage


class TaxNode:  # pylint: disable=too-few-public-methods
    """An NCBI Taxonomy node linked to its parent and children."""

    def __init__(self, tax_id, parent_id):
        self.tax_id = tax_id
        self.parent_id = parent_id
        self.names = []

        self.parent_node = None
        self.families = []
        self.children = []
        self.ancestral = 0


class ClassificationNode:  # pylint: disable=too-few-public-methods
    """A Dfam Classification node linked to its parent and children."""

    def __init__(
        self, class_id, parent_id, name, type_name, subtype_name
    ):  # pylint: disable=too-many-arguments
        self.class_id = class_id
        self.parent_id = parent_id
        self.name = name
        self.type_name = type_name
        self.subtype_name = subtype_name

        self.parent_node = None
        self.children = []

    def full_name(self):
        """
        Returns the full name of this classification node, with the name of each
        classification level delimited by a semicolon.
        """
        name = self.name
        node = self.parent_node

        while node is not None:
            name = node.name + ";" + name
            node = node.parent_node

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
                raise TypeError(
                    "Incompatible type for '{}'. Expected '{}', got '{}'".format(
                        name, expected_type, type(value)
                    )
                ) from exc
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
        return "%s.%s '%s': %s len=%d" % (
            self.accession,
            self.version,
            self.name,
            self.classification,
            self.length or -1,
        )

    def to_dfam_hmm(
        self,
        famdb,
        species=None,
        include_class_in_name=False,
        require_general_threshold=False,
    ):  # pylint: disable=too-many-locals,too-many-branches
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
                try:
                    (hmm_ga, hmm_tc, hmm_nc, hmm_fdr) = map(float, parts[1:])
                except Exception as err:
                    hmm_ga = 0.0
                    hmm_tc = 0.0
                    hmm_nc = 0.0
                    hmm_fdr = 0.0
                    print(
                        "Error in thresholds for accession={} and taxid={}".format(
                            self.accession_with_optional_version(), tax_id
                        ),
                        file=sys.stderr,
                    )

                # only recover name, do need for partition number
                tax_name = famdb.get_taxon_name(tax_id, "scientific name")[0]
                if tax_id == species:
                    species_hmm_ga, species_hmm_tc, species_hmm_nc = (
                        hmm_ga,
                        hmm_tc,
                        hmm_nc,
                    )
                th_lines += [
                    "TaxId:%d; TaxName:%s; GA:%.2f; TC:%.2f; NC:%.2f; fdr:%.3f;"
                    % (tax_id, tax_name, hmm_ga, hmm_tc, hmm_nc, hmm_fdr)
                ]

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
        out += "\n".join(model_lines[i + 1 :])

        return out

    __COMPLEMENT_TABLE = str.maketrans("ACGTRYWSKMNXBDHV", "TGCAYRSWMKNXVHDB")

    def to_fasta(
        self,
        famdb,
        use_accession=False,
        include_class_in_name=False,
        do_reverse_complement=False,
        buffer=None,
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

            sequence = sequence[buffer[0] - 1 : buffer[1]]
            identifier = f"{identifier}#buffer"

        if do_reverse_complement:
            sequence = sequence.translate(self.__COMPLEMENT_TABLE)
            sequence = sequence[::-1]

        if include_class_in_name and not buffer:
            rm_class = self.repeat_type
            if self.repeat_subtype:
                rm_class += "/" + self.repeat_subtype
            identifier = f"{identifier}#{rm_class}"

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
            out += sequence[i : i + 60] + "\n"
            i += 60

        return out

    def to_embl(
        self, famdb, include_meta=True, include_seq=True
    ):  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
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
        append("AC", self.accession + ";")
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
                append(
                    "KW",
                    "Long terminal repeat of retrovirus-like element; %s." % self.name,
                )
            else:
                append(
                    "KW", "%s/%s." % (self.repeat_type or "", self.repeat_subtype or "")
                )
            out += "XX\n"

            for clade_id in self.clades:
                lineage = famdb.get_lineage_path(clade_id, partition=False)
                if lineage[0] == ["root"]:
                    lineage = lineage[1:]
                if len(lineage) > 0:
                    append("OS", lineage[-1])
                    append("OC", "; ".join(lineage[:-1]) + ".", True)
            out += "XX\n"

            if self.citations:
                citations = json.loads(self.citations)
                citations.sort(key=lambda c: c["order_added"])
                for cit in citations:
                    append(
                        "RN", "[%d] (bases 1 to %d)" % (cit["order_added"], self.length)
                    )
                    append("RA", cit["authors"], True)
                    append("RT", cit["title"], True)
                    append("RL", cit["journal"])
                    out += "XX\n"

            append("CC", self.description, True)
            out += "CC\n"
            append("CC", "RepeatMasker Annotations:")
            append("CC", "     Type: %s" % (self.repeat_type or ""))
            append("CC", "     SubType: %s" % (self.repeat_subtype or ""))

            species_names = [famdb.get_sanitized_name(c) for c in self.clades]
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

                    append(
                        "FT",
                        "CDS             %d..%d" % (cds["cds_start"], cds["cds_end"]),
                    )
                    append_featuredata('/product="%s"' % cds["product"])
                    append_featuredata("/number=%s" % cds["exon_count"])
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

            append(
                "SQ",
                "Sequence %d BP; %d A; %d C; %d G; %d T; %d other;"
                % (
                    len(sequence),
                    counts["a"],
                    counts["c"],
                    counts["g"],
                    counts["t"],
                    counts["other"],
                ),
            )

            while i < len(sequence):
                chunk = sequence[i : i + 60]
                i += 60

                j = 0
                line = ""
                while j < len(chunk):
                    line += chunk[j : j + 10] + " "
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
                match = re.match(r"(\S*)", value)
                acc = match.group(1)
                acc = acc.rstrip(";")
                family.accession = acc
            elif code == "NM":
                family.name = value
            elif code == "DE":
                family.description = value
            elif code == "CC":
                # TODO: Consider only recognizing these after seeing "RepeatMasker Annotations"

                matches = re.match(r"\s*Type:\s*(\S+)", value)
                if matches:
                    family.repeat_type = matches.group(1).strip()

                matches = re.match(r"\s*SubType:\s*(\S+)", value)
                if matches:
                    family.repeat_subtype = matches.group(1).strip()

                matches = re.search(r"Species:\s*(.+)", value)
                if matches:
                    for spec in matches.group(1).split(","):
                        name = spec.strip()
                        if name:
                            tax_id = lookup.get(name.lower())
                            if tax_id is not None:
                                family.clades += [tax_id]
                            else:
                                LOGGER.warning("Could not find taxon for '%s'", name)
                matches = re.search(r"SearchStages:\s*(\S+)", value)
                if matches:
                    family.search_stages = matches.group(1).strip()

                matches = re.search(r"BufferStages:\s*(\S+)", value)
                if matches:
                    family.buffer_stages = matches.group(1).strip()

                matches = re.search(r"Refineable", value)
                if matches:
                    family.refineable = True

        header = ""
        family = None
        in_header = True
        in_metadata = False

        nodes = lookup.values()

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
                            header_line = matches.group(2).rstrip("*").strip()
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
                        keep = False
                        for clade in family.clades:
                            if clade in nodes:
                                LOGGER.debug(
                                    f"Including {family.accession} in taxa {clade} from {filename}"
                                )
                                keep = True
                        if keep:
                            yield family
                        family = None

                    # Part of the sequence area
                    else:
                        family.consensus += re.sub(r"[^A-Za-z]", "", line)

        # if header_cb:
        #     header_cb(header)
