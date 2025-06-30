#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Usage: famdb.py [-h] [-l LOG_LEVEL] [-i DB_DIR] command ...

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
    Anthony Gray <anthony.gray@systemsbiology.org>
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
import json
import logging
import os
import re
import sys
import traceback

from famdb_globals import (
    LOGGER,
    FILE_DESCRIPTION,
    FAMILY_FORMATS_EPILOG,
    MISSING_FILE,
    HELP_URL,
)
from famdb_classes import FamDB


# Command-line utilities
def command_info(args):
    """The 'info' command displays some of the stored metadata."""
    db_info = args.db_dir.get_metadata()
    counts = args.db_dir.get_counts()
    print()
    print(
        f"""\
FamDB Directory     : {os.path.realpath(args.db_dir.db_dir)}
FamDB Format Version: {db_info["famdb_version"]}
FamDB Creation Date : {db_info["created"]}

Database: {db_info["name"]}
Version : {db_info["db_version"]}
Date    : {db_info["date"]}

{db_info["description"]}

{counts['file']} Partitions Present
Total consensus sequences present: {counts["consensus"]}
Total HMMs present               : {counts["hmm"]}
"""
    )
    args.db_dir.show_files()
    if args.history:
        args.db_dir.show_history()


def command_names(args):
    """The 'names' command displays all names of all taxa that match the search term."""

    entries = []
    entries += args.db_dir.resolve_names(args.term)

    if args.format == "pretty":
        prev_exact = None
        for tax_id, is_exact, partition, names in entries:
            if is_exact != prev_exact:
                if is_exact:
                    print("Exact Matches\n=============")
                else:
                    if prev_exact:
                        print()
                    print("Non-exact Matches\n=================")
                prev_exact = is_exact

            print(
                f"Taxon: {tax_id}, Partition: {partition}, Names: {', '.join([f'{n[1]} ({n[0]})' for n in names])}"
            )

    elif args.format == "json":
        obj = []
        for tax_id, is_exact, partition, names in entries:
            names_obj = [{"kind": name[0], "value": name[1]} for name in names]
            obj += [{"id": tax_id, "partition": partition, "names": names_obj}]
        print(json.dumps(obj))
    else:
        raise ValueError("Unimplemented names format: %s" % args.format)


def print_lineage_tree(
    file,
    tree,
    gutter_self,
    gutter_children,
    curated_only=False,
    uncurated_only=False,
):
    """Pretty-prints a lineage tree with box drawing characters."""

    if not tree:
        return
    if type(tree) == str:
        tax_id = tree
        children = []
    else:
        tax_id = tree[0]
        children = tree[1:]
    name, tax_partition = file.get_taxon_name(tax_id, "scientific name")
    if name != "Not Found":
        fams = file.get_families_for_taxon(
            tax_id,
            tax_partition,
            curated_only=curated_only,
            uncurated_only=uncurated_only,
        )
        num_fams = len(fams) if fams is not None else 0
        missing_message = MISSING_FILE % (tax_partition, file.db_dir, HELP_URL)
        missing_message = (
            missing_message.replace("\t", f"{gutter_self[:-2]}│ * \t")
            + f"\n{gutter_self[:-2]}│"
        )
        count = f"[{num_fams}]" if fams is not None else missing_message
        print(f"{gutter_self}{tax_id} {name}({tax_partition}) {count}")

    # All but the last child need a downward-pointing line that will link up
    # to the next child, so this is split into two cases
    if len(children) > 1:
        for child in children[:-1]:
            print_lineage_tree(
                file,
                child,
                gutter_children + "├─",
                gutter_children + "│ ",
                curated_only,
                uncurated_only,
            )

    if children:
        print_lineage_tree(
            file,
            children[-1],
            gutter_children + "└─",
            gutter_children + "  ",
            curated_only,
            uncurated_only,
        )


def print_lineage_semicolons(
    file,
    tree,
    parent_name,
    starting_at,
    curated_only=False,
    uncurated_only=False,
):
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
    name, tax_partition = file.get_taxon_name(tax_id, "scientific name")

    if name != "Not Found":
        if parent_name:
            name = parent_name + ";" + name

        if starting_at == tax_id:
            starting_at = None

        if not starting_at:
            fams = file.get_families_for_taxon(
                tax_id, tax_partition, curated_only, uncurated_only
            )
            count = (
                f"[{len(fams)}]"
                if fams is not None
                else f"(Taxon in Partition {tax_partition}, Partition File Not Found)"
            )
            print(f"{tax_id}({tax_partition}): {name} {count}")

        for child in children:
            print_lineage_semicolons(
                file,
                child,
                name,
                starting_at,
                curated_only,
                uncurated_only,
            )


def get_lineage_totals(
    file,
    tree,
    target_id,
    partition,
    curated_only=False,
    uncurated_only=False,
    seen=None,
    present=None,
):
    """
    Recursively calculates the total number of families
    on ancestors and descendants of 'target_id' in the given 'tree'.

    'seen' is required to track families that are present on multiple
    lineages due to horizontal transfer and ensure each family
    is only counted one time, either as an ancestor or a descendant.
    """
    if not seen:
        seen = set()
    if not present:
        present = set()

    tax_id = tree[0]
    children = tree[1:]
    partition = file.find_taxon(tax_id)
    accessions = file.get_families_for_taxon(
        tax_id, partition, curated_only, uncurated_only
    )

    count_here = 0
    for acc in accessions:
        if acc not in seen:
            seen.add(acc)
            count_here += 1

    if target_id == tax_id:
        target_id = None

    counts = [0, 0]
    for child in children:
        partition = file.find_taxon(tax_id)
        if partition is not None:
            new_counts, new_present = get_lineage_totals(
                file,
                child,
                target_id,
                partition,
                curated_only,
                uncurated_only,
                seen,
                present,
            )
            counts[0] += new_counts[0]
            counts[1] += new_counts[1]
            present.add(partition)
            present.update(new_present)

    if target_id is None:
        counts[1] += count_here
    else:
        counts[0] += count_here
    return counts, present


def command_lineage(args):
    """The 'lineage' command outputs ancestors and/or descendants of the given taxon."""

    target_id, partition = args.db_dir.resolve_one_species(args.term)

    if not target_id:
        print(f"No species found for search term '{args.term}'", file=sys.stderr)
        return
    if target_id == "Ambiguous":
        return
    tree = args.db_dir.get_lineage(
        target_id,
        descendants=args.descendants,
        ancestors=args.ancestors or args.format == "semicolon",
        complete=args.complete or args.format == "semicolon",
    )
    if not tree:
        return
    if args.format == "pretty":
        print_lineage_tree(
            args.db_dir,
            tree,
            "",
            "",
            args.curated,
            args.uncurated,
        )
    elif args.format == "semicolon":
        print_lineage_semicolons(
            args.db_dir, tree, "", target_id, args.curated, args.uncurated
        )
    elif args.format == "totals":
        totals, present = get_lineage_totals(
            args.db_dir, tree, target_id, partition, args.curated, args.uncurated
        )
        present = (
            ", ".join([str(val) for val in present]) + ";" if present else partition
        )
        missing = (
            " absent related partitions: "
            + ", ".join([str(val) for val in set(tree.missing.values())])
            if hasattr(tree, "missing")
            else ""
        )
        print(
            f"{totals[0]} entries in ancestors; {totals[1]} lineage-specific entries; found in partitions: {present}{missing}"
        )
    else:
        raise ValueError("Unimplemented lineage format: %s" % args.format)


def print_families(args, families, header, species=None):
    """
    Prints each family in 'families', optionally with a copyright header. The
    format is determined by 'args.format' and additional data (such as
    taxonomy) is taken from 'args.db_dir'.

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
        db_info = args.db_dir.get_metadata()
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
            if include_class_in_name:
                name = family.name or family.accession
                rm_class = family.repeat_type
                if family.repeat_subtype:
                    rm_class += "/" + family.repeat_subtype
                family.name = name + "#" + rm_class
            entry = str(family) + "\n"
        elif args.format == "hmm":
            entry = family.to_dfam_hmm(
                args.db_dir,
                include_class_in_name=include_class_in_name,
                require_general_threshold=require_general_threshold,
            )
        elif args.format == "hmm_species":
            entry = family.to_dfam_hmm(
                args.db_dir,
                species,
                include_class_in_name=include_class_in_name,
                require_general_threshold=require_general_threshold,
            )
        elif (
            args.format == "fasta"
            or args.format == "fasta_name"
            or args.format == "fasta_acc"
        ):
            use_accession = args.format == "fasta_acc"

            buffers = []
            if stage and family.buffer_stages:
                for spec in family.buffer_stages.split(","):
                    if "[" in spec:
                        matches = re.match(r"(\d+)\[(\d+)-(\d+)\]", spec.strip())
                        if matches:
                            if stage == int(matches.group(1)):
                                buffers += [
                                    [int(matches.group(2)), int(matches.group(3))]
                                ]
                        else:
                            LOGGER.warning(
                                "Ingored invalid buffer specification: '%s'",
                                spec.strip(),
                            )
                    else:
                        buffers += [stage == int(spec)]

            if not buffers:
                buffers += [None]

            entry = ""
            for buffer_spec in buffers:
                entry += (
                    family.to_fasta(
                        args.db_dir,
                        use_accession=use_accession,
                        include_class_in_name=include_class_in_name,
                        buffer=buffer_spec,
                    )
                    or ""
                )

                if add_reverse_complement:
                    entry += (
                        family.to_fasta(
                            args.db_dir,
                            use_accession=use_accession,
                            include_class_in_name=include_class_in_name,
                            do_reverse_complement=True,
                            buffer=buffer_spec,
                        )
                        or ""
                    )
        elif args.format == "embl":
            entry = family.to_embl(args.db_dir)
        elif args.format == "embl_meta":
            entry = family.to_embl(args.db_dir, include_meta=True, include_seq=False)
        elif args.format == "embl_seq":
            entry = family.to_embl(args.db_dir, include_meta=False, include_seq=True)
        else:
            raise ValueError("Unimplemented family format: %s" % args.format)

        if entry:
            print(entry, end="")


def command_family(args):
    """The 'family' command outputs a single family by name or accession."""
    family = args.db_dir.get_family_by_accession(args.accession)
    if not family:
        family = args.db_dir.get_family_by_name(args.accession)

    if family:
        print_families(args, [family], False)


def command_families(args):
    """The 'families' command outputs all families associated with the given taxon."""
    target_id, _ = args.db_dir.resolve_one_species(args.term)
    if not target_id:
        print(f"No species found for search term '{args.term}'", file=sys.stderr)
        return
    elif target_id == "Ambiguous":
        return

    families = []

    is_hmm = args.format.startswith("hmm")

    # NB: This is speed-inefficient, because get_accessions_filtered needs to
    # read the whole family data even though we read it again right after.
    # However it is *much* more memory-efficient than loading all the family
    # data at once and then sorting by accession.
    accessions = sorted(
        args.db_dir.get_accessions_filtered(
            tax_id=target_id,
            descendants=args.descendants,
            ancestors=args.ancestors,
            curated_only=args.curated,
            uncurated_only=args.uncurated,
            is_hmm=is_hmm,
            stage=args.stage,
            repeat_type=args.repeat_type,
            name=args.name,
        )
    )
    families = map(args.db_dir.get_family_by_accession, accessions)

    header = True if accessions else False
    print_families(args, families, header, target_id)


# RepeatMasker Commands -----------------------------------------------------------------------
def command_fasta_all(args):
    """
    command prints out all curated families in FASTA format
    This command is not documented in the help. It is used to export all of the curated families
    to FASTA format for use by RepeatMasker
    """
    args.format = "fasta_name"
    args.include_class_in_name = True
    print_families(args, args.db_dir.fasta_all("/DF"), True, 1)
    print_families(args, args.db_dir.fasta_all("/Aux"), True, 1)


def command_repeatpeps(args):
    """prints the RepeatPeps file"""
    print(args.db_dir.get_repeatpeps())


def command_edit_description(args):
    """Updates the db description"""
    args.db_dir.update_description(args.new)


def command_append(args):
    """
    The 'append' command reads an EMBL file and appends its entries to an
    existing famdb file.
    """

    lookup = args.db_dir.get_all_taxa_names()
    # infile_lookup = {}
    # with open(args.infile) as file:
    #     infile_lookup = json.load(file)
    # lookup.update(infile_lookup)

    header = None

    def set_header(val):
        nonlocal header
        header = val

    embl_iter = FamDB.read_embl_families(args.infile, lookup, header_cb=set_header)

    message = f"Adding Families From {args.infile.split('/')[-1]}"
    rec = args.db_dir.append_start_changelog(message)

    LOGGER.info(message)
    total_ctr = 0
    added_ctr = 0
    file_counts = {}
    new_val_taxa = set()
    dups = set()
    missing_files = {}

    for entry in embl_iter:
        # check installation namespace and skip entry if it already exists
        if entry.accession in args.rb_names or not args.db_dir.check_unique(entry):
            LOGGER.info(
                f"Skipped {entry.accession}. A family with the same accession/name is already in Dfam"
            )
            continue

        total_ctr += 1
        acc = entry.accession
        added = False

        # prepare set of local files to add family to
        add_files = set()
        add_taxa = set()
        for clade in entry.clades:
            file = args.db_dir.find_taxon(clade)
            if args.db_dir.files.get(file):
                if args.db_dir.files[file].has_taxon(clade):
                    add_files.add(file)
                    # check if the taxon is empty
                    if not args.db_dir.get_families_for_taxon(clade, file):
                        add_taxa.add(clade)
            else:
                missing_files[file] = missing_files.get(file, 0) + 1

        if not add_files:
            LOGGER.debug(f" {acc} not added to local files, local file not found")

        for file in add_files:
            try:
                args.db_dir.files[file].add_family(entry)
                LOGGER.debug(f"Added {acc} to file {file}")
                if not added:
                    added_ctr += 1
                    added = True
                file_counts[file] = file_counts.get(file, 0) + 1
            except Exception as e:
                LOGGER.debug(f" Ignoring duplicate entry {entry.accession}: {e}")
                dups.add(entry.accession)

        # track formerly empty clades with new additions
        if added:
            new_val_taxa.update(add_taxa)

    args.db_dir.append_finish_changelog(message, rec)
    args.db_dir.update_changelog(added_ctr, total_ctr, file_counts, args.infile)

    LOGGER.info(f"Added {added_ctr}/{total_ctr} families")
    if dups:
        LOGGER.debug(f" {len(dups)} Duplicate Accesisons: {dups}")
    if missing_files:
        for file in missing_files:
            LOGGER.info(
                f"Partition File {file} Not Found. {missing_files[file]} Entries Were Not Appended:"
            )

    db_info = args.db_dir.get_metadata()

    if args.name:
        db_info["name"] = args.name
    if args.description:
        db_info["description"] += "\n" + args.description

    if header:
        db_info["copyright"] += f"\n\n{header}"

    args.db_dir.set_db_info(
        db_info["name"],
        db_info["db_version"],
        db_info["date"],
        db_info["description"],
        db_info["copyright"],
    )

    # Write the updated counts and metadata
    if new_val_taxa:
        LOGGER.info("Rebuilding Sparse Taxonomy Tree")
        args.db_dir.rebuild_pruned_tree(new_val_taxa)

    LOGGER.info("Finalizing Files")
    args.db_dir.finalize()

def build_args():
    """builds and parses the command line args"""
    parser = argparse.ArgumentParser(
        description=FILE_DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-l", "--log_level", default="INFO")

    parser.add_argument("-i", "--db_dir", help="specifies the directory to query")

    subparsers = parser.add_subparsers(
        description="""Specifies the kind of query to perform.
For more information on all the possible options for a command, add the --help option after it:
famdb.py families --help
""",
        #  metavar, if specified overrides what shows up on the help line as valid
        #  subcommands.  All subcommands will however be printed in the error message
        #  if a bad subcommand is entered as a possibility, so it doesn't hide it
        #  completely.  This is added to hide the new fasta_all command.
        metavar="{info,names,lineage,families,family,append}",
    )
    # INFO --------------------------------------------------------------------------------------------------------------------------------
    p_info = subparsers.add_parser(
        "info", description="List general information about the file."
    )
    p_info.add_argument(
        "--history",
        action="store_true",
        help="List the file changelog in addition to general information",
    )
    p_info.set_defaults(func=command_info)

    # NAMES --------------------------------------------------------------------------------------------------------------------------------
    p_names = subparsers.add_parser(
        "names", description="List the names and taxonomy identifiers of a clade."
    )
    p_names.add_argument(
        "-f",
        "--format",
        default="pretty",
        choices=["pretty", "json"],
        metavar="<format>",
        help="choose output format. The default is 'pretty'. 'json' is more appropriate for scripts.",
    )
    p_names.add_argument(
        "term",
        nargs="+",
        help="search term. Can be an NCBI taxonomy identifier or part of a scientific or common name",
    )
    p_names.set_defaults(func=command_names)

    # LINEAGE --------------------------------------------------------------------------------------------------------------------------------
    p_lineage = subparsers.add_parser(
        "lineage",
        description="List the taxonomy tree including counts of families at each clade.",
    )
    p_lineage.add_argument(
        "-a",
        "--ancestors",
        action="store_true",
        help="include all ancestors of the given clade",
    )
    p_lineage.add_argument(
        "-d",
        "--descendants",
        action="store_true",
        help="include all descendants of the given clade",
    )
    p_lineage.add_argument(
        "-k",
        "--complete",
        action="store_true",
        help="include output of taxa without families",
        default=False,
    )
    p_lineage.add_argument(
        "-c",
        "--curated",
        action="store_true",
        help="only tabulate curated families ('DF' records)",
    )
    p_lineage.add_argument(
        "-u",
        "--uncurated",
        action="store_true",
        help="only tabulate uncurated families ('DR' records)",
    )
    p_lineage.add_argument(
        "-f",
        "--format",
        default="pretty",
        choices=["pretty", "semicolon", "totals"],
        metavar="<format>",
        help="choose output format. The default is 'pretty'. 'semicolon' is more appropriate for scripts. 'totals' displays the number of ancestral and lineage-specific families found.",
    )
    p_lineage.add_argument(
        "term",
        nargs="+",
        help="search term. Can be an NCBI taxonomy identifier or an unambiguous scientific or common name",
    )
    p_lineage.set_defaults(func=command_lineage)

    # FAMILIES --------------------------------------------------------------------------------------------------------------------------------
    family_formats = [
        "summary",
        "hmm",
        "hmm_species",
        "fasta_name",
        "fasta_acc",
        "embl",
        "embl_meta",
        "embl_seq",
    ]
    family_formats_epilog = FAMILY_FORMATS_EPILOG

    p_families = subparsers.add_parser(
        "families",
        description="Retrieve the families associated \
with a given clade, optionally filtered by additional criteria",
        epilog=family_formats_epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_families.add_argument(
        "-a",
        "--ancestors",
        action="store_true",
        help="include all ancestors of the given clade",
    )
    p_families.add_argument(
        "-d",
        "--descendants",
        action="store_true",
        help="include all descendants of the given clade",
    )
    p_families.add_argument(
        "--stage",
        type=int,
        help="include only families that should be searched in the given stage",
    )
    p_families.add_argument(
        "--class",
        dest="repeat_type",
        type=str,
        help="include only families that have the specified repeat Type/SubType",
    )
    p_families.add_argument(
        "--name",
        type=str,
        help="include only families whose name begins with this search term",
    )
    p_families.add_argument(
        "-u",
        "--uncurated",
        action="store_true",
        help="include only 'uncurated' families (i.e. named DRXXXXXXXXX)",
    )
    p_families.add_argument(
        "-c",
        "--curated",
        action="store_true",
        help="include only 'curated' families (i.e. not named DFXXXXXXXXX)",
    )
    p_families.add_argument(
        "-f",
        "--format",
        default="summary",
        choices=family_formats,
        metavar="<format>",
        help="choose output format.",
    )
    p_families.add_argument(
        "--add-reverse-complement",
        action="store_true",
        help="include a reverse-complemented copy of each matching family; only suppported for fasta formats",
    )
    p_families.add_argument(
        "--include-class-in-name",
        action="store_true",
        help="include the RepeatMasker type/subtype after the name (e.g. HERV16#LTR/ERVL); only supported for hmm and fasta formats",
    )
    p_families.add_argument(
        "--require-general-threshold",
        action="store_true",
        help="skip families missing general thresholds (and log their accessions at the debug log level)",
    )
    p_families.add_argument(
        "term",
        nargs="+",
        help="search term. Can be an NCBI taxonomy identifier or an unambiguous scientific or common name",
    )
    p_families.set_defaults(func=command_families)

    # FAMILY --------------------------------------------------------------------------------------------------------------------------------
    p_family = subparsers.add_parser(
        "family",
        description="Retrieve details of a single family.",
        epilog=family_formats_epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p_family.add_argument(
        "-f",
        "--format",
        default="summary",
        choices=family_formats,
        metavar="<format>",
        help="choose output format.",
    )
    p_family.add_argument(
        "accession", help="the accession of the family to be retrieved"
    )
    p_family.set_defaults(func=command_family)

    # APPEND --------------------------------------------------------------------------------------------------------------------------------
    p_append = subparsers.add_parser("append")
    p_append.add_argument("infile", help="the name of the input file to be appended")
    p_append.add_argument(
        "exclusion_file",
        help="the name of the file listing family names to be excluded",
    )
    p_append.add_argument(
        "--name", help="new name for the database (replaces the existing name)"
    )
    p_append.add_argument(
        "--description",
        help="additional database description (added to the existing description)",
    )
    p_append.set_defaults(func=command_append)

    # FASTA ALL --------------------------------------------------------------------------------------------------------------------------------
    p_fasta = subparsers.add_parser("fasta_all")
    p_fasta.set_defaults(func=command_fasta_all)

    # RepeatPeps -------------------------------------------------------------------------------------------------------------------------------
    p_rp = subparsers.add_parser("repeat_peps")
    p_rp.set_defaults(func=command_repeatpeps)

    # Edit Description -------------------------------------------------------------------------------------------------------------------------------
    p_desc = subparsers.add_parser("edit_description")
    p_desc.add_argument("new")
    p_desc.set_defaults(func=command_edit_description)

    return parser


def main():  # ================================================================================================================================
    """Parses command-line arguments and runs the requested command."""

    logging.basicConfig()

    parser = build_args()
    args = parser.parse_args()
    logging.getLogger().setLevel(getattr(logging, args.log_level.upper()))

    write_commands = [command_append, command_edit_description]
    if "func" in args and args.func in write_commands:
        mode = "r+"
    else:
        mode = "r"

    if "term" in args:
        args.term = " ".join(args.term)

    # For RepeatMasker: Try Libraries/RepeatMaskerLib.h5, if no file was specified
    # in the arguments and that file exists.
    if not args.db_dir:
        # sys.path[0], if non-empty, is initially set to the directory of the
        # originally-invoked script.
        if sys.path[0]:
            default_db_dir = os.path.join(sys.path[0], "Libraries/famdb")
            if os.path.exists(default_db_dir):
                args.db_dir = default_db_dir

    if not (args.db_dir and os.path.exists(args.db_dir) and os.path.isdir(args.db_dir)):
        LOGGER.error(
            "Please specify a directory containing FamDB files to operate on with the -i/--file option."
        )
        exit(1)

    if hasattr(args,'func') and args.func.__name__ == "command_append":
        if os.path.exists(args.exclusion_file):
            try:
                with open(args.exclusion_file) as f:
                    args.rb_names = set(name.strip() for name in f.readlines())
            except Exception:
                LOGGER.error(f"{args.exclusion_file} could not be parsed.")
            exit(1)
        else:
            LOGGER.error(f"{args.exclusion_file} not found.")
            exit(1)

    try:
        args.db_dir = FamDB(args.db_dir, mode)
    except:
        args.db_dir = None
        raise

    if not args.db_dir:
        return

    if "func" in args:
        try:
            args.func(args)
        except Exception:
            traceback.print_exc()
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
