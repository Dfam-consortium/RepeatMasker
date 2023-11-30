import datetime
import time
import os
import json
import sys

import h5py
import numpy

from famdb_helper_classes import Family, Lineage
from famdb_globals import (
    LOGGER,
    FILE_VERSION,
    GENERATOR_VERSION,
    LEAF_LINK,
    ROOT_LINK,
    GROUP_FAMILIES,
    GROUP_LOOKUP_BYNAME,
    GROUP_LOOKUP_BYSTAGE,
    GROUP_NODES,
    GROUP_TAXANAMES,
    MISSING_FILE,
    HELP_URL,
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
                "Invalid file mode. Expected 'r' or 'r+' or 'w', got '{}'".format(mode)
            )

        self.filename = filename
        self.file = h5py.File(filename, mode)
        self.mode = mode

        try:
            if reading and self.file.attrs["version"] != FILE_VERSION:
                raise Exception(
                    "File version is {}, but this is version {}".format(
                        self.file.attrs["version"],
                        FILE_VERSION,
                    )
                )
        except:
            # This 'except' catches both "version" missing from attrs, or the
            # value not matching if it is present.
            raise Exception("This file cannot be read by this version of famdb.py.")

        if self.mode == "w":
            self.seen = {}
            self.added = {"consensus": 0, "hmm": 0}
            self.__write_metadata()
        elif self.mode == "r+":
            self.added = self.get_counts()

    # Export Setters ----------------------------------------------------------------------------------------------------
    def set_partition_info(self, partition_num):
        """Sets partition number (key to file info) and bool if is root file or not"""
        self.file.attrs["partition_num"] = partition_num
        self.file.attrs["root"] = partition_num == "0" or partition_num == 0

    def set_file_info(self, map_str):
        """Stores information about other files as json string"""
        self.file.attrs["file_info"] = json.dumps(map_str)

    def set_db_info(self, name, version, date, desc, copyright_text):
        """Sets database metadata for the current file"""
        self.file.attrs["db_name"] = name
        self.file.attrs["db_version"] = version
        self.file.attrs["db_date"] = date
        self.file.attrs["db_description"] = desc
        self.file.attrs["db_copyright"] = copyright_text

    def __write_metadata(self):
        """Sets file data during writing"""
        self.file.attrs["generator"] = f"famdb.py v{GENERATOR_VERSION}"
        self.file.attrs["version"] = FILE_VERSION
        self.file.attrs["created"] = str(datetime.datetime.now())

    def finalize(self):
        """Writes some collected metadata, such as counts, to the database"""
        self.file.attrs["count_consensus"] = self.added["consensus"]
        self.file.attrs["count_hmm"] = self.added["hmm"]

    # Attribute Getters -----------------------------------------------------------------------------------------------
    def get_partition_num(self):
        """Partition num is used as the key in file_info"""
        return self.file.attrs["partition_num"]

    def get_file_info(self):
        """returns dictionary containing information regarding other related files"""
        return json.loads(self.file.attrs["file_info"])

    def is_root(self):
        """Tests if file is root file"""
        return self.file.attrs["root"]

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
        'generator', 'version', 'created', 'partition_name', 'partition_detail'
        """
        num = self.file.attrs["partition_num"]
        partition = self.get_file_info()["file_map"][str(num)]
        return {
            "generator": self.file.attrs["generator"],
            "version": self.file.attrs["version"],
            "created": self.file.attrs["created"],
            "partition_name": partition["T_root_name"],
            "partition_detail": ", ".join(partition["F_roots_names"]),
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

    # File Utils
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

        # TODO: This is awkward. The EMBL files being appended may only have an
        # "accession", but that accession may match the *name* of a family
        # already in Dfam. The accession may also match a family already in
        # Dfam, but with a "v" added.

        # check by accession first
        accession = family.accession
        binned_acc = accession_bin(accession)
        binned_v = accession_bin(accession + 'v')

        if self.file.get(f"{binned_acc}/{accession}") or self.file.get(f"{binned_v}/{accession}v"):
            return False

        # check for unique name
        #if family.name:
        #    name_lookup = f"{GROUP_LOOKUP_BYNAME}/{family.name}"
        #    if self.file.get(name_lookup) or self.file.get(name_lookup + 'v'):
        #        return False

        if self.file.get(f"{GROUP_LOOKUP_BYNAME}/{accession}") or self.file.get(f"{GROUP_LOOKUP_BYNAME}/{accession}v"):
            return False

        return True


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
            self.file.require_group(GROUP_LOOKUP_BYNAME)[
                str(family.name)
            ] = h5py.SoftLink(fam_link)
        # In FamDB format version 0.5 we removed the /Families/ByAccession group as it's redundant
        # (all the data is in Families/<datasets> *and* HDF5 suffers from poor performance when
        # the number of entries in a group exceeds 200-500k.

        for clade_id in family.clades:
            clade = str(clade_id)
            nodes = self.file[GROUP_NODES]
            if clade in nodes:
                families_group = nodes[clade].require_group("Families")
                families_group[family.accession] = h5py.SoftLink(fam_link)

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

        LOGGER.debug("Added family %s (%s)", family.name, family.accession)

    # Taxonomy Nodes
    def write_taxonomy(self, tax_db, nodes):
        """Writes taxonomy nodes in 'nodes' to the database."""
        LOGGER.info("Writing taxonomy nodes")
        start = time.perf_counter()

        count = 0
        for node in nodes:
            count += 1
            group = self.file.require_group(GROUP_NODES).require_group(
                str(tax_db[node].tax_id)
            )
            parent_id = int(tax_db[node].parent_id) if node != 1 else None
            if parent_id:
                group.create_dataset("Parent", data=numpy.array([parent_id]))

            child_ids = []
            for child in tax_db[node].children:
                child_ids += [int(child.tax_id)]
            group.create_dataset("Children", data=numpy.array(child_ids))
        delta = time.perf_counter() - start
        LOGGER.info("Wrote %d taxonomy nodes in %f", count, delta)

    # Data Access Methods ------------------------------------------------------------------------------------------------
    def has_taxon(self, tax_id):
        """Returns True if 'self' has a taxonomy entry for 'tax_id'"""
        return str(tax_id) in self.file[GROUP_NODES]

    def get_families_for_taxon(self, tax_id, curated_only=False, uncurated_only=False):
        """Returns a list of the accessions for each family directly associated with 'tax_id'."""
        group = (
            self.file[GROUP_NODES][str(tax_id)].get("Families")
            if f"{GROUP_NODES}/{tax_id}/Families" in self.file
            else {}
        )

        # Filter out DF/DR or not at all depending on flags
        if curated_only:
           return list(filter(lambda x: (x[1] == 'F'), group.keys()))
        elif uncurated_only:
            return list(filter(lambda x: (x[1] == 'R'), group.keys()))
        else:
            return list(group.keys())


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
        root = self.is_root()
        if descendants:

            def descendants_of(tax_id):
                descendants = [int(tax_id)]
                for child in group_nodes[str(tax_id)]["Children"]:
                    # only list the decendants of the target node if it's not being combined with another decendant lineage
                    if not kwargs.get("for_combine") and str(child) in group_nodes:
                        descendants += [descendants_of(child)]
                    elif root:
                        descendants += [f"{LEAF_LINK}{child}"]
                return descendants

            tree = descendants_of(tax_id)
        else:
            tree = [tax_id]

        if ancestors:
            while tax_id:
                node = group_nodes[str(tax_id)]
                if "Parent" in node:
                    if str(node["Parent"][0]) in group_nodes:
                        tax_id = node["Parent"][0]
                        tree = [tax_id, tree]
                    else:
                        tree = [f"{ROOT_LINK}{tax_id}", tree]
                        tax_id = None
                else:
                    tax_id = None

        lineage = Lineage(tree, root, self.get_partition_num())
        return lineage

    def filter_stages(self, accession, stages):
        """Returns True if the family belongs to a search or buffer stage in 'stages'."""
        for stage in stages:
            grp = self.file[GROUP_LOOKUP_BYSTAGE].get(stage)
            if grp and accession in grp:
                return True

        return False

    # Family Getters --------------------------------------------------------------------------
    def get_family_names(self):  # TODO unused
        """Returns a list of names of families in the database."""
        return sorted(self.file[GROUP_LOOKUP_BYNAME].keys(), key=str.lower)

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
        entry = self.file[GROUP_LOOKUP_BYNAME].get(name)
        return get_family(entry)


class FamDBRoot(FamDBLeaf):
    def __init__(self, filename, mode="r"):
        super(FamDBRoot, self).__init__(filename, mode)

        if mode == "r" or mode == "r+":
            self.names_dump = {
                partition: json.loads(
                    self.file[f"{GROUP_TAXANAMES}/{partition}"]["TaxaNames"][0]
                )
                for partition in self.file[GROUP_TAXANAMES]
            }
            self.file_info = self.get_file_info()
            self.__lineage_cache = {}

    def write_taxa_names(self, tax_db, nodes):
        """
        Writes Names -> taxa maps per partition
        """
        LOGGER.info("Writing TaxaNames")
        for partition in nodes:
            taxnames_group = self.file.require_group(GROUP_TAXANAMES + f"/{partition}")
            names_dump = {}
            for node in nodes[partition]:
                names_dump[node] = tax_db[node].names
            names_data = numpy.array([json.dumps(names_dump)])
            names_dset = taxnames_group.create_dataset(
                "TaxaNames", shape=names_data.shape, dtype=FamDBLeaf.dtype_str
            )
            names_dset[:] = names_data

    def get_taxon_names(self, tax_id):
        """
        Checks names_dump for each partition and returns a list of [name_class, name_value, partition]
        of the taxon given by 'tax_id'.
        """
        for partition in self.names_dump:
            names = self.names_dump[partition].get(str(tax_id))
            if names:
                return names
        return []

    def get_taxon_name(self, tax_id, kind="scientific name"):
        """
        Checks names_dump for each partition and returns eturns the first name of the given 'kind'
        for the taxon given by 'tax_id', or None if no such name was found.
        """
        for partition in self.names_dump:
            names = self.names_dump[partition].get(str(tax_id))
            if names is not None:
                for name in names:
                    if name[0] == kind:
                        return [name[1], int(partition)]
        return "Not Found", "N/A"

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
        for partition in self.names_dump:
            for tax_id, names in self.names_dump[partition].items():
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
                    yield [int(tax_id), exact, int(partition)]

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
            for partition in self.names_dump:
                if str(tax_id) in self.names_dump[partition]:
                    return [[tax_id, int(partition), True]]

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
                for tax_id, _ in similar_results:
                    names = self.get_taxon_names(tax_id)
                    print(
                        tax_id,
                        ", ".join(["{1}".format(*n) for n in names]),
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
        """

        name = self.get_taxon_name(tax_id, "scientific name")
        if name:
            name = sanitize_name(name[0])
        return name

    def get_lineage_path(self, tax_id, tree=[], cache=True, partition=True):
        """
        Returns a list of strings encoding the lineage for 'tax_id'.
        """

        if cache and tax_id in self.__lineage_cache:
            return self.__lineage_cache[tax_id]
        if not tree:
            tree = self.get_lineage(tax_id, ancestors=True)

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
        for partition in self.names_dump:
            if str(tax_id) in self.names_dump[partition]:
                return int(partition)
        return None

    def parent_of(self, tax_id):
        group_nodes = self.file[GROUP_NODES]
        for node in group_nodes:
            if int(tax_id) in group_nodes[node]["Children"]:
                return node
        return None


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
        db_prefixes = {}
        h5_files = []
        for file in os.listdir(db_dir):
            if file.endswith(".h5"):
                h5_files += [file]
            if file.endswith(".0.h5"):
                db_prefixes[file[:-5]] = 1

        # Make sure we only have at least one database present
        if len(db_prefixes) == 0:
            if h5_files:
                LOGGER.error("A partitioned famdb database is not present in " + db_dir + "\n" + \
                             "There were several *.h5 files present however, they do not appear\n" + \
                             "to be in the correct format: " + "\n".join(h5_files) + "\n")
            else:
                LOGGER.error("A partitioned famdb database is not present in " + db_dir )
            exit(1)

        # Make sure we have *only* one database present
        if len(db_prefixes) > 1:
            LOGGER.error("Multiple famdb root partitions were found in this export directory: " + \
                          ", ".join(db_prefixes.keys()) + "\nEach famdb database " + \
                          "should be in separate folders.")
            exit(1)

        # Tabulate all partitions for db_prefix
        db_prefix = list(db_prefixes.keys())[0]
        for file in h5_files:
            if db_prefix in file:
                fields = file.split(".")
                idx = int( fields[-2] )
                if idx == 0:
                    self.files[idx] = FamDBRoot(f"{db_dir}/{file}", mode)
                else:
                    self.files[idx] = FamDBLeaf(f"{db_dir}/{file}", mode)

        file_info = self.files[0].get_file_info()

        self.db_dir = db_dir
        self.file_map = file_info["file_map"]
        self.uuid = file_info["meta"]["partition_id"]
        self.db_version = file_info["meta"]["db_version"]
        self.db_date = file_info["meta"]["db_date"]

        err_files = []
        for file in self.files:
            meta = self.files[file].get_file_info()["meta"]
            if (
                self.uuid != meta["partition_id"]
                or self.db_version != meta["db_version"]
                or self.db_date != meta["db_date"]
            ):
                err_files += [file]
        if err_files:
            LOGGER.error(f"Files From Different Partitioning Runs: {err_files}")
            exit()

    def get_lineage_combined(self, tax_id, **kwargs):
        # check if tax_id exists in Dfam
        location = self.find_taxon(tax_id)
        if location is None:
            print("Taxon Not Found In Dfam")
            return None
        if location not in self.files:
            print(MISSING_FILE % (location, self.db_dir, HELP_URL))
            return None
        # query lineage in correct file
        base_lineage = self.files[location].get_lineage(tax_id, **kwargs)

        if base_lineage.descendants:  # lineage extends from root file to leaf file(s)
            add_lineages = []
            missing = {}
            for taxa in base_lineage.links[LEAF_LINK].values():
                # find location of each linked node
                loc = self.find_taxon(taxa)
                if loc and loc in self.files:
                    # query and save subtree if file is installed
                    add_lineages += [
                        self.files[loc].get_lineage(
                            taxa, descendants=True, ancestors=True
                        )
                    ]
                elif loc and loc not in self.files:
                    # if file is not found, return Lineage of 1 node, record that node is missing
                    add_lineages += [
                        Lineage([f"{ROOT_LINK}{taxa}", [int(taxa)]], False, loc)
                    ]
                    missing[taxa] = loc

            # Combine lineages
            for lin in add_lineages:
                base_lineage += lin

            # attach missing file info
            base_lineage.missing = missing

        if base_lineage.ancestors:  # lineage extends from leaf file to root file
            # find ancestor node in root and query lineage
            ancestor_node = self.files[0].parent_of(
                list(base_lineage.links[ROOT_LINK].keys())[0]
            )  # TODO this is probably really slow
            root_lineage = self.files[0].get_lineage(
                ancestor_node, descendants=True, ancestors=True, for_combine=True
            )
            base_lineage += root_lineage

        # strip out leftover links
        def remove_links(lineage):
            for thing in list(lineage):
                if not thing or type(thing) == str:
                    lineage.remove(thing)
                if type(thing) == list:
                    remove_links(thing)

        if kwargs.get("remove_links") or base_lineage.descendants:
            remove_links(base_lineage)

        return base_lineage

    def show_files(self):
        # repbase_file = "./partitions/RMRB_spec_to_tax.json" TODO
        print(f"\nPartition Details\n-----------------")
        for part in sorted([int(x) for x in self.file_map]):
            part_str = str(part)
            partition_name = self.file_map[part_str]["T_root_name"]
            partition_detail = ', '.join(self.file_map[part_str]["F_roots_names"])
            filename = self.file_map[part_str]["filename"]
            if part in self.files:
                print(f" Partition {part} [{filename}]: {partition_name} {f'- {partition_detail}' if partition_detail else ''}")
                counts = self.files[part].get_counts()
                print(f"     Consensi: {counts['consensus']}, HMMs: {counts['hmm']}")
            else:
                print(f" Partition {part} [ Absent ]: {partition_name} {f'- {partition_detail}' if partition_detail else ''}")
            print()


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
                # yield from self.file[FamDBLeaf.GROUP_LOOKUP_BYACC].keys() # TODO unused
                for file in files:
                    names = families_iterator(
                        files[file].file[GROUP_FAMILIES], "Families"
                    )
                    for name in names:
                        yield name
            else:
                lineage = self.get_lineage_combined(
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

    def resolve_names(self, term):
        entries = []
        for tax_id, partition, is_exact in self.files[0].resolve_species(term):
            names = self.files[0].get_taxon_names(tax_id)
            entries += [[tax_id, is_exact, partition, names]]
        return entries

    def fasta_all(self, group):
        seen = set()
        for file in self.files:
            if GROUP_FAMILIES+group in self.files[file].file:
                for name in families_iterator(self.files[file].file[GROUP_FAMILIES+group], "Families" + group):
                    if name not in seen:
                        seen.add(name)
                        yield self.get_family_by_accession(name)

    # Wrapper methods ---------------------------------------------------------------------------------------
    def get_counts(self):
        counts = {"consensus": 0, "hmm": 0, "file": 0}
        for file in self.files:
            file_counts = self.files[file].get_counts()
            counts["consensus"] += file_counts["consensus"]
            counts["hmm"] += file_counts["hmm"]
            counts["file"] += 1
        return counts

    def get_lineage_path(self, tax_id, **kwargs):
        """method used in EMBL exports"""
        lineage = self.get_lineage_combined(tax_id, **kwargs)
        partition = (
            kwargs.get("partition") if kwargs.get("partition") is not None else True
        )
        cache = kwargs.get("cache") if kwargs.get("cache") is not None else True
        return self.files[0].get_lineage_path(
            tax_id, lineage, cache=cache, partition=partition
        )

    def get_sanitized_name(self, tax_id):
        """method used in EMBL exports"""
        return self.files[0].get_sanitized_name(tax_id)

    def get_db_info(self):
        return self.files[0].get_db_info()

    def resolve_one_species(self, term):
        return self.files[0].resolve_one_species(term)

    def get_metadata(self):
        return self.files[0].get_metadata()

    def get_taxon_name(self, tax_id, kind):
        return self.files[0].get_taxon_name(tax_id, kind)

    def get_families_for_taxon(self, tax_id, partition, curated_only=False,
                               uncurated_only=False):
        if partition in self.files:
            return self.files[partition].get_families_for_taxon(tax_id, curated_only, uncurated_only)
        else:
            return None

    def get_family_by_accession(self, accession):
        for file in self.files:
            fam = self.files[file].get_family_by_accession(accession)
            if fam:
                return fam
        return None

    def get_family_by_name(self, accession):
        for file in self.files:
            fam = self.files[file].get_family_by_name(accession)
            if fam:
                return fam
        return None

    def find_taxon(self, tax_id):
        return self.files[0].find_taxon(tax_id)

    def finalize(self):
        for file in self.files:
            self.files[file].finalize()

    def set_db_info(self, name, version, date, desc, copyright_text):
        for file in self.files:
            self.files[file].set_db_info(name, version, date, desc, copyright_text)

    def filter_stages(self, accession, stages):
        for file in self.files:
            fam = self.files[file].get_family_by_accession(accession)
            if fam:
                return self.files[file].filter_stages(accession, stages)

    # File Utils
    def close(self):
        """Closes this FamDB instance, making further use invalid."""
        for file in self.files:
            self.files[file].close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
