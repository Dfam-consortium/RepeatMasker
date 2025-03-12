import datetime
import time
import os
import json
import sys
import re
import h5py
import numpy

from famdb_helper_classes import Family, TaxNode
from famdb_globals import (
    LOGGER,
    FAMDB_VERSION,
    GROUP_FAMILIES,
    GROUP_LOOKUP_BYNAME,
    GROUP_LOOKUP_BYSTAGE,
    GROUP_LOOKUP_BYTAXON,
    GROUP_NODES,
    GROUP_FILE_HISTORY,
    GROUP_REPEATPEPS,
    DATA_CHILDREN,
    DATA_PARENT,
    DATA_VAL_CHILDREN,
    DATA_VAL_PARENT,
    DATA_TAXANAMES,
    DATA_PARTITION,
    DATA_NAMES_CACHE,
    META_DB_VERSION,
    META_DB_DESCRIPTION,
    META_DB_COPYRIGHT,
    META_DB_DATE,
    META_DB_NAME,
    META_CREATED,
    META_META,
    META_UUID,
    META_FILE_INFO,
    META_FAMDB_VERSION,
    META_FILE_MAP,
    DESCRIPTION,
)
from famdb_helper_methods import (
    sanitize_name,
    sounds_like,
    families_iterator,
    filter_curated,
    filter_repeat_type,
    filter_search_stages,
    filter_name,
    get_family,
    accession_bin,
    is_fasta,
)


class FamDBLeaf:
    """Transposable Element Family and taxonomy database."""

    dtype_str = h5py.special_dtype(vlen=str)

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
            raise ValueError(
                f"Invalid file mode. Expected 'r' or 'r+' or 'w', got '{mode}'"
            )

        self.filename = filename
        self.file = h5py.File(filename, mode)
        self.mode = mode

        if (reading and not self.file.attrs.get(META_FAMDB_VERSION)) or (
            reading and not self.version_match()
        ):
            LOGGER.error(
                f"\t Partition {self.get_partition_num()}:This file cannot be read by this version of famdb.py.\n"
                f" Export File Version: {self.file.attrs.get(META_FAMDB_VERSION, 'Not Found')}\n"
                f" FamDB Script Version: {FAMDB_VERSION}\n"
            )
            sys.exit(1)

        if self.mode == "w":
            self.seen = {}
            self.added = {"consensus": 0, "hmm": 0}
            self.__write_metadata()
        elif self.mode == "r+":
            self.added = self.get_counts()

    def version_match(self):
        file_version = self.file.attrs.get(META_FAMDB_VERSION)
        file_splits = file_version.split(".")
        file_major = file_splits[0] if file_splits else None

        script_splits = FAMDB_VERSION.split(".")
        script_major = script_splits[0] if script_splits else None

        same_major = file_major == script_major

        if not same_major:
            return False
        return True

    def update_changelog(self, message, verified=False):
        """
        Creates a OtherData/FileHistory/Timestamp/Message/bool
        to record file changes. Defaults to False to show that change is not complete
        """
        time_stamp = str(datetime.datetime.now())
        group = self.file.require_group(GROUP_FILE_HISTORY).require_group(time_stamp)
        group.create_dataset(message, data=numpy.array([verified]))
        return time_stamp

    def _verify_change(self, time_stamp, message):
        """
        Sets the data of a log entry to True, indicating that it was successful
        """
        self.file[GROUP_FILE_HISTORY][time_stamp][message][0] = True

    def _change_logger(func):
        """
        A wrapper method to update and verify the changelog for common methods
        """
        func_to_note = {
            "__write_metadata": "File Initialized",
            "set_metadata": "Metadata Set",
            "add_family": "Family Added",
            "write_repeatpeps": "RepeatPeps Written",
            "write_taxonomy": "Taxonomy Nodes Written",
            "write_full_taxonomy": "Taxonomy Nodes Written",
            "update_description": "File Description Updated",
            "update_pruned_taxa": "Pruned Tree Updated",
        }
        message = func_to_note[func.__name__]

        def wrapper(self, *args, **kwargs):
            time_stamp = self.update_changelog(message)
            func(self, *args, **kwargs)
            self._verify_change(time_stamp, message)

        return wrapper

    # Export Setters ----------------------------------------------------------------------------------------------------
    @_change_logger
    def __write_metadata(self):
        """Sets file data during writing. Called during file creation"""
        self.file.attrs[META_FAMDB_VERSION] = FAMDB_VERSION
        self.file.attrs[META_CREATED] = str(datetime.datetime.now())
        self.file.attrs[META_DB_DESCRIPTION] = DESCRIPTION

    @_change_logger
    def set_metadata(self, partition_num, map_str, name, version, date, copyright_text):
        """
        Sets database metadata for the current file
        Stores information about other files as json string
        Sets partition number (key to file info) and bool if is root file or not
        """
        self.file.attrs[META_DB_NAME] = name
        self.file.attrs[META_DB_VERSION] = version
        self.file.attrs[META_DB_DATE] = date
        self.file.attrs[META_DB_COPYRIGHT] = copyright_text

        self.file.attrs[META_FILE_INFO] = json.dumps(map_str)

        self.file.attrs["partition_num"] = partition_num
        self.file.attrs["root"] = partition_num == "0" or partition_num == 0

    def finalize(self):
        """Writes some collected metadata, such as counts, to the database"""
        self.file.attrs["count_consensus"] = self.added["consensus"]
        self.file.attrs["count_hmm"] = self.added["hmm"]

    @_change_logger
    def update_description(self, new_desc):
        """Updates the description. Available to the user and during the append command"""
        self.file.attrs[META_DB_DESCRIPTION] = new_desc

    # Attribute Getters -----------------------------------------------------------------------------------------------
    def get_partition_num(self):
        """Partition num is used as the key in file_info"""
        return self.file.attrs["partition_num"]

    def get_file_info(self):
        """returns dictionary containing information regarding other related files"""
        return json.loads(self.file.attrs[META_FILE_INFO])

    def is_root(self):
        """Tests if file is root file"""
        return self.file.attrs["root"]

    def get_metadata(self):
        """
        Gets file metadata for the current file as a dict with keys
        'famdb_version', 'created', 'partition_name', 'partition_detail',
        'db_name', 'db_version', 'db_date', 'db_description', 'db_copyright'
        """
        if "db_name" not in self.file.attrs:
            return None
        num = self.get_partition_num()
        partition = self.get_file_info()[META_FILE_MAP][str(num)]
        return {
            "famdb_version": self.file.attrs[META_FAMDB_VERSION],
            "created": self.file.attrs[META_CREATED],
            "partition_name": partition["T_root_name"],
            "partition_detail": ", ".join(partition["F_roots_names"]),
            "name": self.file.attrs[META_DB_NAME],
            "db_version": self.file.attrs[META_DB_VERSION],
            "date": self.file.attrs[META_DB_DATE],
            "description": self.file.attrs[META_DB_DESCRIPTION],
            "copyright": self.file.attrs[META_DB_COPYRIGHT],
        }

    def get_history(self):
        """
        Retrieves and concatenates the changelog into a string
        """
        history = self.file.get(GROUP_FILE_HISTORY)
        messages = {stamp: list(history[stamp].keys())[0] for stamp in history.keys()}
        hist_str = f"\n File {self.get_partition_num()}\n"
        for entry in messages:
            hist_str += f"{entry} - {messages[entry]}\n"
        return hist_str

    def get_counts(self):
        """
        Gets counts of entries in the current file as a dict
        with 'consensus', 'hmm'
        """
        return {
            "consensus": self.file.attrs["count_consensus"],
            "hmm": self.file.attrs["count_hmm"],
        }

    # File Utils
    def interrupt_check(self):
        """
        Changelogs Start as False and are flipped to True when complete
        Returns bool if any changes are not confirmed
        """
        interrupted = False
        history = self.file.get(GROUP_FILE_HISTORY)
        for el in history:
            item = history.get(el)
            note = list(item.keys())[0]
            val = item[note][()][0]
            if not val:
                interrupted = True
                break
        return interrupted

    def close(self):
        """Closes this FamDB instance, making further use invalid."""
        self.file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    # Data Writing Methods ---------------------------------------------------------------------------------------------
    # Family Methods
    def __check_unique(self, family):
        """Verifies that 'family' is uniquely identified by its value of 'key'."""

        # This is awkward. The EMBL files being appended may only have an
        # "accession", but that accession may match the *name* of a family
        # already in Dfam. The accession may also match a family already in
        # Dfam, but with a "v" added.
        # This has been spot-checked and seems to avoid conflicts - Anthony 11/5/24

        # check by accession first
        accession = family.accession
        binned_acc = accession_bin(accession)
        binned_v = accession_bin(accession + "v")

        if self.file.get(f"{binned_acc}/{accession}") or self.file.get(
            f"{binned_v}/{accession}v"
        ):
            return False

        # check for unique name
        # if family.name:
        #    name_lookup = f"{GROUP_LOOKUP_BYNAME}/{family.name}"
        #    if self.file.get(name_lookup) or self.file.get(name_lookup + 'v'):
        #        return False

        if self.file.get(f"{GROUP_LOOKUP_BYNAME}/{accession}") or self.file.get(
            f"{GROUP_LOOKUP_BYNAME}/{accession}v"
        ):
            return False

        return True

    # no @_change_logger here to avoid 1000s of history logs. it is called in the methods that call add_family
    def add_family(self, family):
        """Adds the family described by 'family' to the database."""
        # Verify uniqueness of name and accession.
        # This is important because of the links created to them later.
        if not self.__check_unique(family):
            raise Exception(
                f"Family is not unique! Already seen {family.accession} {f'({family.name})' if family.name else ''}"
            )

        # Increment counts
        if family.consensus:
            self.added["consensus"] += 1
        if family.model:
            self.added["hmm"] += 1

        # Create the family data
        # In v0.5 we bin the datasets into subgroups to improve performance
        group_path = accession_bin(family.accession)
        dset = self.file.require_group(group_path).create_dataset(
            family.accession, (0,)
        )

        # Set the family attributes
        for k in Family.META_LOOKUP:
            value = getattr(family, k)
            if value:
                dset.attrs[k] = value

        # Create links
        fam_link = f"/{group_path}/{family.accession}"
        if family.name:
            self.file.require_group(GROUP_LOOKUP_BYNAME)[str(family.name)] = (
                h5py.SoftLink(fam_link)
            )
        # In FamDB format version 0.5 we removed the /Families/ByAccession group as it's redundant
        # (all the data is in Families/<datasets> *and* HDF5 suffers from poor performance when
        # the number of entries in a group exceeds 200-500k.

        for clade_id in family.clades:
            clade = str(clade_id)
            nodes = self.file[GROUP_LOOKUP_BYTAXON]
            if clade in nodes:
                nodes[clade][family.accession] = h5py.SoftLink(fam_link)

        def add_stage_link(stage, accession):
            stage_group = self.file.require_group(GROUP_LOOKUP_BYSTAGE).require_group(
                stage.strip()
            )
            if accession not in stage_group:
                stage_group[accession] = h5py.SoftLink(fam_link)

        if family.search_stages:
            for stage in family.search_stages.split(","):
                add_stage_link(stage, family.accession)

        if family.buffer_stages:
            for stage in family.buffer_stages.split(","):
                stage = stage.split("[")[0]
                add_stage_link(stage, family.accession)

        LOGGER.debug(f"Added family {family.name} ({family.accession})")

    # Taxonomy Nodes
    @_change_logger
    def write_taxonomy(self, nodes):
        """Writes taxonomy nodes to the database. These nodes only contain links to family data and not data regarding tree relationships"""
        LOGGER.info(f"Writing taxonomy in partition")
        start = time.perf_counter()

        count = 0
        for node in nodes:
            count += 1
            self.file.require_group(GROUP_LOOKUP_BYTAXON).require_group(str(node))
        delta = time.perf_counter() - start
        LOGGER.info(f"Wrote {count} taxonomy nodes in {delta}")

    # Data Access Methods ------------------------------------------------------------------------------------------------
    def has_taxon(self, tax_id):
        """Returns True if 'self' has a taxonomy entry for 'tax_id'"""
        return str(tax_id) in self.file[GROUP_LOOKUP_BYTAXON]

    def get_families_for_taxon(self, tax_id, curated_only=False, uncurated_only=False):
        """Returns a list of the accessions for each family directly associated with 'tax_id'."""
        group = (
            self.file[GROUP_LOOKUP_BYTAXON][str(tax_id)]
            if f"{GROUP_LOOKUP_BYTAXON}/{tax_id}" in self.file
            else {}
        )

        # Filter out DF/DR or not at all depending on flags
        if curated_only:
            return list(filter(lambda a: filter_curated(a, True), group.keys()))
        elif uncurated_only:
            return list(filter(lambda a: filter_curated(a, False), group.keys()))
        else:
            return list(group.keys())

    def filter_stages(self, accession, stages):
        """Returns True if the family belongs to a search or buffer stage in 'stages'."""
        for stage in stages:
            grp = self.file[GROUP_LOOKUP_BYSTAGE].get(stage)
            if grp and accession in grp:
                return True

        return False

    # Family Getters --------------------------------------------------------------------------
    def get_family_by_accession(self, accession):
        """Returns the family with the given accession."""
        path = accession_bin(accession)
        if path in self.file:
            entry = self.file[path].get(accession)
            return get_family(entry)
        return None

    def get_family_by_name(self, name):
        """Returns the family with the given name."""
        # TODO: This will also suffer the performance issues seen with
        #       other groups that exceed 200-500k entries in a single group
        #       at some point.  This needs to be refactored to scale appropriately.
        # There are 24,768 names as of Dfam 3.8 - Anthony
        entry = self.file[GROUP_LOOKUP_BYNAME].get(name)
        return get_family(entry)


class FamDBRoot(FamDBLeaf):
    def __init__(self, filename, mode="r"):
        super(FamDBRoot, self).__init__(filename, mode)

        if mode == "r" or mode == "r+":
            self.names_dump = json.loads(self.file[DATA_NAMES_CACHE][()].decode())
            self.file_info = self.get_file_info()
            self.__lineage_cache = {}

    @FamDBLeaf._change_logger
    def write_full_taxonomy(self, tax_db, nodes):
        """
        Takes a map of TaxaNodes
        Writes taxonomy nodes to the database.
        Includes parent-child relationships
        Also cache all taxa names as a node:[names] json string
        This cache is loaded on __init__ to speed up search times
        """
        LOGGER.info(f"Writing Full Taxonomy Tree Root File")
        start = time.perf_counter()

        partition_map = {
            node: int(partition) for partition in nodes for node in nodes[partition]
        }
        names_dump = {}
        count = 0
        for node in tax_db:
            count += 1
            group = self.file.require_group(GROUP_NODES).require_group(
                str(tax_db[node].tax_id)
            )
            parent_id = int(tax_db[node].parent_id) if node != 1 else None
            if parent_id:
                group.create_dataset(DATA_PARENT, data=numpy.array([parent_id]))

            child_ids = []
            for child in tax_db[node].children:
                child_ids += [int(child.tax_id)]
            group.create_dataset(DATA_CHILDREN, data=numpy.array(child_ids))

            names = tax_db[node].names
            group.create_dataset(DATA_TAXANAMES, data=numpy.array(names, dtype="S"))
            names_dump[node] = names
            group.create_dataset(
                DATA_PARTITION, data=numpy.array([partition_map[node]])
            )

        LOGGER.info(f"Writing Name Cache String")
        self.file.create_dataset(
            DATA_NAMES_CACHE, data=numpy.array(json.dumps(names_dump), dtype="S")
        )

        delta = time.perf_counter() - start
        LOGGER.info(f"Wrote {count} taxonomy nodes in full tree in {delta}")

    @FamDBLeaf._change_logger
    def update_pruned_taxa(self, tree):
        """
        Takes a map of TaxaNodes
        Updates the nodes to include sparse parent-child relationships
        based on which nodes have family data associated with them
        """
        for id in tree:
            node = tree[id]
            val_children = [int(child) for child in node.val_children]
            val_parent = int(node.val_parent) if node.val_parent else None
            group = self.file[GROUP_NODES][str(id)]
            if group.get(DATA_VAL_CHILDREN):
                del group[DATA_VAL_CHILDREN]
            group.create_dataset(
                DATA_VAL_CHILDREN,
                data=numpy.array(val_children),
                shape=(len(val_children),),
                dtype="i8",
            )
            if val_parent:
                if group.get(DATA_VAL_PARENT):
                    del group[DATA_VAL_PARENT]
                group.create_dataset(
                    DATA_VAL_PARENT,
                    data=numpy.array([val_parent]),
                    shape=(1,),
                    dtype="i8",
                )

    @FamDBLeaf._change_logger
    def write_repeatpeps(self, infile):
        """
        Writing RepeatPeps to its own group as one big string.
        For now, only RepeatModeler consumes this, and does so
        by loading the whole file, so no need to do more
        """
        LOGGER.info(f"Writing RepeatPeps From File: {infile}")
        fasta = is_fasta(infile)
        if fasta:
            with open(infile, "r") as file:
                repeatpeps_str = file.read()
                rp_data = self.file.create_dataset(
                    GROUP_REPEATPEPS, shape=1, dtype=h5py.string_dtype()
                )
                rp_data[:] = repeatpeps_str
            LOGGER.info("RepeatPeps Saved")
        else:
            LOGGER.error(f"File {infile} not in FASTA format, write cancelled")

    def get_repeatpeps(self):
        """
        Retrieve RepeatPeps File
        """
        return self.file.get(GROUP_REPEATPEPS)[0].decode(
            encoding="UTF-8", errors="strict"
        )

    # currently unused:
    # def get_family_names(self):
    #     """Returns a list of names of families in the database."""
    #     return sorted(self.file[GROUP_LOOKUP_BYNAME].keys(), key=str.lower)

    def get_taxon_names(self, tax_id):
        """
        Checks names_dump for each partition and returns a list of [name_class, name_value, partition]
        of the taxon given by 'tax_id'.
        """
        nodes = self.file[GROUP_NODES]
        node = nodes.get(str(tax_id))
        if node:
            return [
                [name.decode() for name in name_pair]
                for name_pair in node[DATA_TAXANAMES][:]
            ]
        return []

    def get_taxon_name(self, tax_id, kind="scientific name"):
        """
        Checks names_dump for each partition and returns eturns the first name of the given 'kind'
        for the taxon given by 'tax_id', or None if no such name was found.
        """
        failure = ("Not Found", "N/A")

        nodes = self.file[GROUP_NODES]
        node = nodes.get(str(tax_id))
        if not node:
            return failure

        names = [
            [name.decode() for name in name_pair]
            for name_pair in node[DATA_TAXANAMES][:]
        ]
        partition = node[DATA_PARTITION][0]

        if names and partition is not None:
            for name in names:
                if name[0] == kind:
                    return [name[1], partition]
        return failure

    def search_taxon_names(self, text, kind=None, search_similar=False):
        """
        Searches 'self' for taxons with a name containing 'text', returning an
        iterator that yields a tuple of (id, is_exact, partition) for each matching node.
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
                partition = self.find_taxon(tax_id)
                yield [int(tax_id), exact, partition]

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
            if str(tax_id) in self.names_dump:
                partition = self.find_taxon(tax_id)
                return [[tax_id, partition, True]]

            return []
        except ValueError:
            pass

        # Perform a search by name, splitting between exact and inexact matches for sorting
        exact = []
        inexact = []
        for tax_id, is_exact, partition in self.search_taxon_names(
            term, kind, search_similar
        ):
            hit = [tax_id, partition]
            if is_exact:
                exact += [hit]
            else:
                inexact += [hit]

        # Combine back into one list, with exact matches first
        results = [[*hit, True] for hit in exact]
        for hit in inexact:
            results += [[*hit, False]]

        if len(results) == 0 and not search_similar:
            # Try a sounds-like search (currently soundex)
            similar_results = self.resolve_species(term, kind, True)
            if similar_results:
                print(
                    "No results were found for that name, but some names sound similar:",
                    file=sys.stderr,
                )
                for tax_id, partition, exact in similar_results:
                    names = self.get_taxon_names(tax_id)
                    print(
                        tax_id,
                        ", ".join([f"{n}" for n in names]),
                        file=sys.stderr,
                    )

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

        for result in results:  # result -> [tax_id, partition, exact]
            if result[2]:
                exact_matches += [[result[0], result[1]]]
        if len(exact_matches) == 1:
            return exact_matches[0]

        if len(results) == 1:
            return results[0][:2]
        elif len(results) > 1:
            print(
                f"""Ambiguous search term '{term}' (found {len(results)} results, {len(exact_matches)} exact).
Please use a more specific name or taxa ID, which can be looked
up with the 'names' command.""",
                file=sys.stderr,
            )
            return "Ambiguous", "Ambiguous"
        return None, None

    def get_sanitized_name(self, tax_id):
        """
        Returns the "sanitized name" of tax_id, which is the sanitized version
        of the scientific name.
        Used in EMBL exports
        """

        name = self.get_taxon_name(tax_id, "scientific name")
        if name:
            name = sanitize_name(name[0])
        return name

    def get_lineage(self, tax_id, **kwargs):
        """
        Returns the lineage of 'tax_id'. Recognized kwargs: 'descendants' to include
        descendant taxa, 'ancestors' to include ancestor taxa.
        IDs are returned as a nested list, for example
        [ 1, [ 2, [3, [4]], [5], [6, [7]] ] ]
        where '2' may have been the passed-in 'tax_id'.
        """

        group_nodes = self.file[GROUP_NODES]
        ancestors = True if kwargs.get("ancestors") else False
        descendants = True if kwargs.get("descendants") else False
        children_key = (
            DATA_VAL_CHILDREN if not kwargs.get("complete") else DATA_CHILDREN
        )
        parent_key = DATA_VAL_PARENT if not kwargs.get("complete") else DATA_PARENT

        if descendants:

            def descendants_of(tax_id):
                descendants = [
                    int(tax_id)
                ]  # h5py is based on numpy, need to cast numpy base64 to python int for serialization in Lineage class
                for child in group_nodes[str(tax_id)][children_key]:
                    descendants += [descendants_of(child)]
                return descendants

            tree = descendants_of(tax_id)
        else:
            tree = [tax_id]

        if ancestors:
            while tax_id:
                node = group_nodes[str(tax_id)]
                if parent_key in node:
                    tax_id = node[parent_key][0]
                    tree = [
                        int(tax_id),
                        tree,
                    ]  # h5py is based on numpy, need to cast numpy base64 to python int for serialization in Lineage class
                else:
                    tax_id = None

        return tree

    def get_lineage_path(self, tax_id, cache=True, partition=True, complete=False):
        """
        Returns a list of strings encoding the lineage for 'tax_id'.
        """

        if cache and tax_id in self.__lineage_cache:
            return self.__lineage_cache[tax_id]
        tree = self.get_lineage(tax_id, ancestors=True, complete=complete)
        lineage = []

        while tree:
            node = tree[0]
            if len(tree) > 1:
                found = False
                for t in tree[1:]:
                    if type(t) == list:
                        tree = t
                        found = True
                        break
                if not found:
                    tree = None
            else:
                tree = None

            tax_name = self.get_taxon_name(node, "scientific name")
            if not partition:
                tax_name = tax_name[0]
            lineage += [tax_name]

        if cache:
            self.__lineage_cache[tax_id] = lineage

        return lineage

    def find_taxon(self, tax_id):
        """
        Returns the partition number containing the taxon
        """
        node = self.file[GROUP_NODES].get(str(tax_id))
        if node:
            return int(node[DATA_PARTITION][0])
        return None

    def get_all_taxa_names(self):
        """
        Returns all taxa names in database.
        Names are cached as taxa : [[name type, name], ...]
        Names are returned as {sanitized_lowercase_name: taxa...}
        Used for mapping EMBL file names to taxa nodes
        Used in append command
        """
        return {
            name[1].lower(): taxon
            for taxon, names in self.names_dump.items()
            for name in names
            if name[0] == "sanitized scientific name" or name[0] == "sanitized synonym"
        }


class FamDB:

    def __init__(self, db_dir, mode):
        """
        Initialize from a directory containing a *partitioned* famdb dataset
        """
        self.files = {}

        ## First, identify if there are any root partitions of a partitioned
        ## famdb in this directory:
        # A partioned famdb file is named *.#.h5 where
        # the number represents the partition number and
        # at a minimum partitition 0 must be present.
        root_prefixes = set()
        prefixes = set()
        h5_files = []
        for file in os.listdir(db_dir):
            if file.endswith(".h5"):
                h5_files += [file]
                splits = file.split(".")
                prefix = ".".join(splits[:-2])
                prefixes.add(prefix)
                if file.endswith(".0.h5"):
                    root_prefixes.add(prefix)

        # ensure only one FamDB export per folder
        if len(prefixes) != 1:
            LOGGER.error("Only one export of FamDB should be present in " + db_dir)
            exit(1)

        # Make sure we only have at least one database present
        if len(root_prefixes) == 0:
            if h5_files:
                LOGGER.error(
                    "A partitioned famdb database is not present in "
                    + db_dir
                    + "\n"
                    + "FamDB requires exactly one root file. There were several *.h5 files present. However, they do not appear\n"
                    + "to be in the correct format: "
                    + "\n".join(h5_files)
                    + "\n"
                )
            else:
                LOGGER.error(
                    "A partitioned famdb root file is not present in " + db_dir
                )
            exit(1)

        # Make sure we have *only* one database present
        if len(root_prefixes) > 1:

            LOGGER.error(
                "Multiple famdb root partitions were found in this export directory: "
                + ", ".join(root_prefixes.keys())
                + "\nEach famdb database "
                + "should be in separate folders."
            )
            exit(1)

        # Tabulate all partitions for db_prefix
        db_prefix = list(root_prefixes)[0]
        for file in h5_files:
            if db_prefix in file:
                fields = file.split(".")
                idx = int(fields[-2])
                if idx == 0:
                    self.files[idx] = FamDBRoot(f"{db_dir}/{file}", mode)
                else:
                    self.files[idx] = FamDBLeaf(f"{db_dir}/{file}", mode)

        file_info = self.files[0].get_file_info()
        self.db_dir = db_dir
        self.file_map = file_info[META_FILE_MAP]
        self.uuid = file_info[META_META][META_UUID]
        self.db_version = file_info[META_META][META_DB_VERSION]
        self.db_date = file_info[META_META][META_DB_DATE]

        partition_err_files = []
        for file in self.files:
            meta = self.files[file].get_file_info()[META_META]
            if (
                self.uuid != meta[META_UUID]
                or self.db_version != meta[META_DB_VERSION]
                or self.db_date != meta[META_DB_DATE]
            ):
                partition_err_files += [file]
        if partition_err_files:
            LOGGER.error(
                f"Files From Different Partitioning Runs: {partition_err_files}"
            )
            exit()

        change_err_files = []
        for file in self.files:
            interrupted = self.files[file].interrupt_check()
            if interrupted:
                change_err_files += [file]
        if change_err_files:
            LOGGER.error(f"Files Interrupted During Edit: {change_err_files}")
            exit()

    # Data writing methods ---------------------------------------------------------------------------------------
    def build_pruned_tree(self):
        """
        Establishes a sparse taxonomy tree where parent-child relationships are restricted to
        nodes with associated family data. For example, a node will be assigned a sparse parent
        as the closest ancestor node with data, rather than its actual parent node, if its
        actual parent node is empty.
        This method exists in FamDB instead of FamDBRoot because it is subject to change after an append,
        and because the associated data is stored in FamDBLeaf files

        Taxonomy Tree is stored as a dictionary of TaxNodes ( self.files[0].file[GROUP_NODES][node] )
            node:
                tax_id: int
                parent_id: int
                val: bool
                children: [int]
                val_parent: int
                val_children: [int]

        If adding a family this should be easily modified by:
             1. Identify the node which the new family is assigned to (or more than one for multiple clades)
             2. Set the val flag to True
             3. For each child in val_children
                 ....
             Go over this with Anthony
        """

        def traverse_val_parents(tree, id):
            """Recurse up the tree ancestor by ancestor until it finds the nearest ancestor with data"""
            node = tree[id]
            if node.parent_id:
                parent = tree.get(node.parent_id)
                if parent:
                    if parent.val:
                        return parent.tax_id
                    else:
                        return traverse_val_parents(tree, parent.tax_id)
            else:
                return None

        def traverse_val_children(tree, id, node_id):
            """
            Adds node to it's parent's list of sparse children
            Continues recursion until it finds an ancestor with data
            """
            node = tree[id]
            if node.parent_id:
                parent = tree.get(node.parent_id)
                if parent:
                    parent.val_children += [node_id]
                    if not parent.val:
                        traverse_val_children(tree, parent.tax_id, node_id)

        LOGGER.info("Reading Taxonomy Tree")
        # read taxonomy tree
        tree = {
            node: self.files[0].file[GROUP_NODES][node]
            for node in self.files[0].file[GROUP_NODES]
        }
        LOGGER.info("Mapping Nodes To Files")
        nodes = {
            file: list(self.files[file].file[GROUP_LOOKUP_BYTAXON].keys())
            for file in [file for file in self.files.keys()]
        }
        LOGGER.info("Determining Which Nodes Have Associated Families")
        # build set of nodes with associated family data
        vals = set(
            [
                id
                for file in nodes
                for id in nodes[file]
                if bool(self.files[file].file[GROUP_LOOKUP_BYTAXON][id].keys())
            ]
        )

        # build TaxNodes in tree
        for id in tree:
            node = tree[id]
            children = node[DATA_CHILDREN][()] if node[DATA_CHILDREN].size > 0 else []
            parent = (
                node[DATA_PARENT][()][0]
                if node.get(DATA_PARENT) and node[DATA_PARENT].size > 0
                else None
            )
            val = id in vals

            tree_node = TaxNode(id, str(parent) if parent else None)
            tree_node.val = val
            tree_node.children = children
            tree[id] = tree_node

        LOGGER.info("Full Tree Prepared")
        # assign each node a val_parent
        for id in tree:
            node = tree[id]
            node.val_parent = traverse_val_parents(tree, node.tax_id)

        # add each node with a value to it's parents as a val_child
        for id in tree:
            node = tree[id]
            if node.val:
                traverse_val_children(tree, node.tax_id, node.tax_id)

        LOGGER.info("Pruned Tree Prepared")

        # update database nodes
        self.files[0].update_pruned_taxa(tree)
        LOGGER.info("Pruned Tree Written")

    def rebuild_pruned_tree(self, new_val_taxa):
        """
        This method takes a list/set of taxon ids that did not have families associated with them,
        but do now due to a recent append command. It resets the val_parent/val_child links in the
        taxonomy tree to account for the fact that there is new data in the tree.
        It assumes that a subject node's val_parent and all ancestor nodes between them will set it
        as one of thier val_children in place of any val_children that it used to share with the
        subject node.
        Likewise, it assumes that any of it's val_children and all child nodes betweem will replace
        thier val_parents with the subject node.
        """

        def build_taxa_node(id, value=False):
            """Builds a TaxNode object from HDF5 data"""
            node = self.files[0].file[GROUP_NODES][str(id)]
            children = node[DATA_CHILDREN][()] if node[DATA_CHILDREN].size > 0 else []
            parent = (
                node[DATA_PARENT][()][0]
                if node.get(DATA_PARENT) and node[DATA_PARENT].size > 0
                else None
            )
            val_children = (
                node[DATA_VAL_CHILDREN][()] if node[DATA_VAL_CHILDREN].size > 0 else []
            )
            val_parent = (
                node[DATA_VAL_PARENT][()][0]
                if node.get(DATA_VAL_PARENT) and node[DATA_VAL_PARENT].size > 0
                else None
            )

            tree_node = TaxNode(id, str(parent) if parent else None)
            tree_node.val = value
            tree_node.children = children
            tree_node.val_children = val_children
            tree_node.val_parent = val_parent

            return tree_node

        # RMH: This parameter default pattern "foo=[]" is dangerous.  The
        #      list generated is global and gets reused between independent
        #      invocations!
        #def climb_non_val_parents(node, ancestor_path=[]):
        #    """collects the nodes between a node and it's val_parent, not inclusive"""
        #    if node.parent_id != node.val_parent:
        #        parent_node = build_taxa_node(node.parent_id)
        #        ancestor_path += [parent_node]
        #        climb_non_val_parents(parent_node, ancestor_path)
        #    return ancestor_path

        def climb_non_val_parents(target_id, node, ancestor_path=None):
            """Collects TaxNodes between a given node and a ancestral
               node defined by target_id (exclusive). """

            if ancestor_path == None:
                ancestor_path = []

            if hasattr(node,"parent_id") and str(node.parent_id) != str(target_id):
                parent_node = build_taxa_node(node.parent_id)
                ancestor_path += [parent_node]
                ancestor_path = climb_non_val_parents(target_id, parent_node, ancestor_path)
            return ancestor_path

        update_nodes = {}
        for id in new_val_taxa:
            node = build_taxa_node(id, value=True)

            seen_tax_ids = set()
            # Fix the descendants of the target node
            for val_child in node.val_children:
                child_node = build_taxa_node(val_child, value=True)

                # Fix the val_child node itself
                if child_node.tax_id not in seen_tax_ids:
                    seen_tax_ids.add(child_node.tax_id)
                    child_node.val_parent = id
                    update_nodes[child_node.tax_id] = child_node

                # Fix the ancestors up until (but not including) the target node
                for ancestor in climb_non_val_parents(id, child_node):
                    if ancestor.tax_id not in seen_tax_ids:
                        seen_tax_ids.add(ancestor.tax_id)
                        ancestor.val_parent = id
                        update_nodes[ancestor.tax_id] = ancestor

            # Gather all nodes above the target node up until its val_parent
            change_ancestors = [build_taxa_node(node.val_parent, value=True)]
            change_ancestors += climb_non_val_parents(node.val_parent, node)

            # All nodes above it should point to it as well, instead of any of its val_children
            for ansc_node in change_ancestors:
                # remove any val_children that are below this node
                for vid in node.val_children:
                    ansc_node.val_children = ansc_node.val_children[ansc_node.val_children != vid]
                # add this node to the ancestral val_children
                ansc_node.val_children = numpy.append(ansc_node.val_children, id)
                update_nodes[ansc_node.tax_id] = ansc_node

            # update the tree for each newly val'd taxon, to avoid tangling pointers when multiple updates occur on the same path
            self.files[0].update_pruned_taxa(update_nodes)
            update_nodes = {}

    def set_db_info(self, name, version, date, desc, copyright_text):
        """Method for resetting metadata"""
        for file in self.files:
            partition_num = self.files[file].get_partition_num()
            file_info = self.files[file].get_file_info()
            self.files[file].set_metadata(
                partition_num,
                file_info,
                name,
                version,
                date,
                copyright_text,
            )
            self.files[file].update_description(desc)

    def append_start_changelog(self, message):
        """
        Called when an append command starts
        """
        rec = {}
        for file in self.files:
            time_stamp = self.files[file].update_changelog(message)
            rec[file] = time_stamp
        return rec

    def append_finish_changelog(self, message, rec):
        """
        Called when an append command finishes successfully
        """
        for file in rec:
            self.files[file]._verify_change(rec[file], message)

    def update_changelog(self, added_ctr, total_ctr, file_counts, infile):
        """Used to add a context log after an append command"""
        filename = infile.split("/")[-1]
        for file in self.files:
            if file in file_counts:
                self.files[file].update_changelog(
                    f"Added {file_counts[file]} of {total_ctr} Families From {filename}",
                    verified=True,
                )
            else:
                self.files[file].update_changelog(
                    f"Found No Relevant Families From {filename}", verified=True
                )
            if file == 0:
                self.files[file].update_changelog(
                    f"Total Families {added_ctr} of {total_ctr} Added To Local Files From {filename}",
                    verified=True,
                )

    # Data access methods ---------------------------------------------------------------------------------------
    def show_files(self):
        """Method to show file information by partition and if those files are present"""
        print(f"\nPartition Details\n-----------------")
        for part in sorted([int(x) for x in self.file_map]):
            part_str = str(part)
            partition_name = self.file_map[part_str]["T_root_name"]
            partition_detail = ", ".join(self.file_map[part_str]["F_roots_names"])
            filename = self.file_map[part_str]["filename"]
            if part in self.files:
                print(
                    f" Partition {part} [{filename}]: {partition_name} {f'- {partition_detail}' if partition_detail else ''}"
                )
                counts = self.files[part].get_counts()
                print(f"     Consensi: {counts['consensus']}, HMMs: {counts['hmm']}")
            else:
                print(
                    f" Partition {part} [ Absent ]: {partition_name} {f'- {partition_detail}' if partition_detail else ''}"
                )
            print()

    def show_history(self):
        """Iterates over all present files and prints each history"""
        print(f"\nFile History\n-----------------")
        for file in self.files:
            print(self.files[file].get_history())

    def get_counts(self):
        """Method gets collected counts from each file present"""
        counts = {"consensus": 0, "hmm": 0, "file": 0}
        for file in self.files:
            file_counts = self.files[file].get_counts()
            counts["consensus"] += file_counts["consensus"]
            counts["hmm"] += file_counts["hmm"]
            counts["file"] += 1
        return counts

    def assemble_filters(self, **kwargs):
        """Define family filters (logically ANDed together)"""
        filters = []
        if kwargs.get("curated_only"):
            filters += [lambda a, f: filter_curated(a, True)]
        if kwargs.get("uncurated_only"):
            filters += [lambda a, f: filter_curated(a, False)]

        filter_stage = kwargs.get("stage")
        stages = []
        if filter_stage:
            if filter_stage == 80:
                # "stage 80" = "all stages", so skip filtering
                pass
            elif filter_stage == 95:
                # "stage 95" = this specific stage list:
                stages = ["35", "50", "55", "60", "65", "70", "75"]
                filters += [lambda a, f: self.filter_stages(a, stages)]
            else:
                stages = [str(filter_stage)]
                filters += [lambda a, f: self.filter_stages(a, stages)]

        # HMM only: add a search stage filter to "un-list" families that were
        # allowed through only because they match in buffer stage
        if kwargs.get("is_hmm") and stages:
            filters += [lambda a, f: filter_search_stages(f(), stages)]

        repeat_type = kwargs.get("repeat_type")
        if repeat_type:
            repeat_type = repeat_type.lower()
            filters += [lambda a, f: filter_repeat_type(f(), repeat_type)]

        name = kwargs.get("name")
        if name:
            name = name.lower()
            filters += [lambda a, f: filter_name(f(), name)]

        return filters, stages, repeat_type, name

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

        filters, stages, repeat_type, name_filter = self.assemble_filters(**kwargs)

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
            files = self.files
            if (
                tax_id == 1
                and descendants
                and stages
                and not repeat_type
                and not name_filter
            ):
                for stage in stages:
                    for file in files:
                        by_stage = files[file].file.get(GROUP_LOOKUP_BYSTAGE)
                        if by_stage:
                            grp = by_stage.get(stage)
                            if grp:
                                yield from grp.keys()

            # special case: Searching the whole database, going directly via
            # Families/ is faster than repeatedly traversing the tree
            elif tax_id == 1 and descendants:
                for file in files:
                    names = families_iterator(
                        files[file].file[GROUP_FAMILIES], GROUP_FAMILIES
                    )
                    for name in names:
                        yield name
            else:
                lineage = self.get_lineage(
                    tax_id, ancestors=ancestors, descendants=descendants
                )
                for node in walk_tree(lineage):
                    location = self.find_taxon(node)
                    fams = self.get_families_for_taxon(node, location)
                    if fams:
                        yield from fams

        for accession in iterate_accs():
            if accession in seen:
                continue
            seen.add(accession)

            cached_family = None

            def family_getter():
                nonlocal cached_family
                if not cached_family:
                    path = accession_bin(accession)
                    for file in self.files:
                        if self.files[file].file.get(path):
                            fam = self.files[file].file[path].get(accession)
                            if fam:
                                cached_family = fam
                return cached_family

            match = True
            for filt in filters:
                if not filt(accession, family_getter):
                    match = False
            if match:
                yield accession

    def fasta_all(self, group):
        """
        Method collects all families in a group
        Used to output all curated data from a db
        """
        seen = set()
        for file in self.files:
            if GROUP_FAMILIES + group in self.files[file].file:
                for name in families_iterator(
                    self.files[file].file[GROUP_FAMILIES + group],
                    GROUP_FAMILIES + group,
                ):
                    if name not in seen:
                        seen.add(name)
                        yield self.get_family_by_accession(name)

    # Root Wrapper methods ---------------------------------------------------------------------------------------
    def resolve_names(self, term):
        """Method to find names matching the search term and map them to the correct file"""
        entries = []
        for tax_id, partition, is_exact in self.files[0].resolve_species(term):
            names = self.files[0].get_taxon_names(tax_id)
            entries += [[tax_id, is_exact, partition, names]]
        return entries

    def get_lineage_path(self, tax_id, **kwargs):
        """method used in EMBL exports"""
        partition = (
            kwargs.get("partition") if kwargs.get("partition") is not None else True
        )
        cache = kwargs.get("cache") if kwargs.get("cache") is not None else True
        complete = (
            kwargs.get("complete") if kwargs.get("complete") is not None else True
        )
        return self.files[0].get_lineage_path(
            tax_id, cache=cache, partition=partition, complete=complete
        )

    def get_sanitized_name(self, tax_id):
        """Wrapper method for the Root get_sanitized_name method"""
        return self.files[0].get_sanitized_name(tax_id)

    def get_lineage(self, tax_id, **kwargs):
        """Wrapper method for the Root get_lineage method"""
        return self.files[0].get_lineage(tax_id, **kwargs)

    def resolve_one_species(self, term):
        """Wrapper method for the Root resolve_one_species method"""
        return self.files[0].resolve_one_species(term)

    def get_metadata(self):
        """Wrapper method for the Root get_metadata method"""
        return self.files[0].get_metadata()

    def get_taxon_name(self, tax_id, kind):
        """Wrapper method for the Root get_taxon_name method"""
        return self.files[0].get_taxon_name(tax_id, kind)

    def find_taxon(self, tax_id):
        """Wrapper method for the Root find_taxon method"""
        return self.files[0].find_taxon(tax_id)

    def get_all_taxa_names(self):
        """Wrapper method for the Root get_all_taxa_names method"""
        return self.files[0].get_all_taxa_names()

    def get_repeatpeps(self):
        """Wrapper method for the Root get_repeatpeps method"""
        return self.files[0].get_repeatpeps()

    # Leaf Wrapper methods ---------------------------------------------------------------------------------------
    def get_families_for_taxon(
        self, tax_id, partition, curated_only=False, uncurated_only=False
    ):
        """Wrapper method to call the Leaf get_families_for_taxon on a specific file"""
        if partition in self.files:
            return self.files[partition].get_families_for_taxon(
                tax_id, curated_only, uncurated_only
            )
        else:
            return None

    def get_family_by_accession(self, accession):
        """Wrapper method to call the Leaf get_family_by_accession"""
        for file in self.files:
            fam = self.files[file].get_family_by_accession(accession)
            if fam:
                return fam
        return None

    def get_family_by_name(self, accession):
        """Wrapper method to call the Leaf get_family_by_name"""
        for file in self.files:
            fam = self.files[file].get_family_by_name(accession)
            if fam:
                return fam
        return None

    def finalize(self):
        """Wrapper method to call the Leaf finalize"""
        for file in self.files:
            self.files[file].finalize()

    def filter_stages(self, accession, stages):
        """Wrapper method to call the Leaf filter_stages"""
        for file in self.files:
            fam = self.files[file].get_family_by_accession(accession)
            if fam:
                return self.files[file].filter_stages(accession, stages)

    def update_description(self, new_desc):
        """Wrapper method to call the Leaf update_description"""
        for file in self.files:
            self.files[file].update_description(new_desc)

    # File Utils
    def close(self):
        """Closes this FamDB instance, making further use invalid."""
        for file in self.files:
            self.files[file].close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    # This method is here because famdb_data_loaders.py imports dfamorm, which is not available to users
    @staticmethod
    def read_embl_families(filename, lookup, header_cb=None):
        """
        This method is here because famdb_data_loaders.py imports dfamorm, which is not available to users

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
                            tax_id = lookup.get(name)
                            if tax_id is not None:
                                family.clades += [tax_id]
                            else:
                                name = name.replace("[", "")
                                name = name.replace("]", "")
                                tax_id = lookup.get(name.lower())
                                if tax_id is not None:
                                    family.clades += [tax_id]
                                else:
                                    LOGGER.warning(
                                        f"Could not find taxon for '{name}' upper or lower: line={value}, and ID={family.accession}"
                                    )
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
