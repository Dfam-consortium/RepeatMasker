#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Usage: ./RM2Bed.py [--help] [--out_dir <path>]
                       [--log_level <number>]
                       [--out_prefix <string>]
                       [--sort_criterion 'family'|'class'|
                           'subclass'|'size'|'diverge']
                       [--split 'family'|'class'|'subclass']
                       [--min_length <number>]
                       [--min_hit_num <number>]
                       [--max_divergence <number>]
                       [--min_divergence <number>]
                       [--ovlp_resolution 'higher_score'|
                          'longer_element'|'lower_divergence']

                       <*.align> or <*.out>

    This Python3 script reads in RepeatMasker a single *.out or
    *.align file and creates one or more BED files with
    optional filtering/grouping.  NOTE: This script requires
    the Pandas python package.

    This script is based on buildSummary.pl from
    the RepeatMasker package, and overlap.py/RM2bed.py
    from David Ray's lab.

    Args:
        --help, -h  : show this help message and exit.
        --out_dir   : directory to store create resulting
                        BED files in.  Default=current
                        working directory.
        --out_prefix: prefix for all output filenames.
                        Default=input file prefix
        --log_level : verbosity of log messages.

    Overlap Resolution:
      RepeatMasker uses a variety of methods to resolve
      overlapping annotation including: information provided
      by flanking elements, age of the family, score
      differentials, and the class of the elements.  There
      are cases where RepeatMasker cannot automatically
      determine which is the best call and does print out
      overlapping annotations ( with the '*' flag at the
      end of the *.out files to indicate overlap has
      occured ).

      This utility can be used to arbitarily resolve
      overlaps in RepeatMasker output.  This is useful
      when it is important to assign every base of the
      input sequence to only one family/class for further
      analysis.  Given the following two types of overlap:

        Basic Overlap:
              |----element1-----------|
                             |---------element2--------|

        Containment:
              |---------element1------------------|
                   |---------element2--------|
           or
              |---------element1------------------|
              |---------element2------------------|

      The overlaps can be resolved in 1 of 5 ways:

        --ovlp_resolution 'higher_score'
            Basic Overlap:
              Overlapping bases are assigned to the element
              with the higher score or the first if they are equal.
            Containtment:
              The lower scoring element is discarded even if it's
              the containing alignment.  If both are equal the
              first one is retained.

        --ovlp_resolution 'longer_element'
            Basic Overlap:
              Overlapping bases are assigned to the longest element
              or first if they are equal.
            Containtment:
              The containing element is retained or the first if they
              are equal.

        --ovlp_resolution 'lower_divergence'
            Basic Overlap:
              Overlapping bases are assigned to the lower divergence
              element.  Currently this uses the Kimura diveregence and
              thus this only works with a *.align file.
            Containment:
              The lower divergence element is retained.  If both
              are equal the first one is retained.
            Simple Repeats are a special case.  They do not have a
            valid Kimura divergence calculated.  The default here is
            to give them a maximum divergence ( e.g the lowest priority ).


      Note that these resolutions are arbitrary and may degrade the
      quality of the overall call for a given locus.  Here we
      have a pair of L2 fragments (joined by RepeatMasker using ID
      numbers) with an unresoloved overlapping alignment on one
      L2 fragment.  If we resolve this in favor of the FalsePostive
      element we are ignoring the evidence provided by the second
      L2 fragment in the vicinity:

      |----L2,ID=1--------|                        |------L2,ID=1---|
                   |-----FalsePositive,ID=2------|

      For large-scale family/class counts this isn't going to have a
      large impact but for fine-scale locus annotation this will not
      always make the correct choices.

      TODO: It may be wiser to use the annotations from the *.out file
      instead of the *.align file.  This is becaus they have been already
      vetted and there is far fewer overalapping elements.  To do this,
      *and* maintain the use of Kimura divergence  one could read the
      *.align file first, average the divergence between alignments with
      the same RepeatMasker "ID" ( the identifier used to link fragments ),
      and store this in a hash using the ID.  Then read and use the *.out
      annotations and append the *.align divergences as needed.  The hard
      part would be to figure out how this might go wrong -- are there
      reasons why you might want to keep the divergence of joined fragments
      different?


SEE ALSO: buildSummary.pl
          RepeatMasker: http://www.repeatmasker.org

AUTHOR(S):
    David Ray
    Robert Hubley <rhubley@systemsbiology.org>

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
#
# Style Guide: See Google Python Style Guide
#     https://github.com/google/styleguide/blob/gh-pages/pyguide.md
#

#
# Module imports
#
import sys
import os
import io
import re
import time
import datetime
import logging
import argparse
import gzip
from operator import itemgetter, attrgetter
import pandas as pd

LOGGER = logging.getLogger(__name__)


def _usage():
    """
    _usage() - Print out docstring for this file

    Args:
        None

    Returns:
        Calls help/pydoc with screen pager on this script file
        and exits normally afterwards.
    """
    # Call to help/pydoc with scriptname ( sans path and file extension )
    help(os.path.splitext(os.path.basename(__file__))[0])
    sys.exit(0)

#
# Other subroutines here
#

def openOptGzipFile(filename, modes='r'):
    """
    openOptGzipFile(filename, modes='r') - Open file with optional gzip
    compression

    Args:
        filename:  The full path to a file
        modes   :  File open modes 'r', 'rw', or 'w'

    Returns:
        A file object
    """
    f = open(filename, 'b'+modes)
    # First two byte signature of a gzip'd file
    if (f.read(2) == b'\x1f\x8b'):
        f.seek(0)
        return io.TextIOWrapper(gzip.GzipFile(fileobj=f), encoding='utf-8')
    else:
        f.seek(0)
        return io.TextIOWrapper(f, encoding='utf-8')


def resolve_using_lower_divergence( cluster ):
    """
    resolve_using_lower_divergence( cluster )

    Resolve overlaps by prioritizing elements
    with a lower kimura divergence.  This
    treats simple repeats in a special way.
    Currently Kimura divergence isn't calculated
    for simple repeats.  This method will
    always prioritize interspersed repeats over
    any simple repeat. [ greedy ]

    e.g. Resolve the following cluster of overlaps
           --------30%-------------
                    -20%-----
                          ----5%----
    as:
                    -20%--
                          ----5%----
    """
    prev_cres = None
    # Sort by divergence ascending.  Treat simple repeats as if they are the
    # highest possible divergence
    cluster = sorted(cluster, key=lambda annot: annot[16] if annot[16] >= 0 else 101)
    for i in range(len(cluster)):
        a_res = cluster[i]
        a_start = a_res[5]
        a_end = a_res[6]
        if ( a_start == 0 and a_end == 0 ):
            continue
        a_len = a_end - a_start + 1
        a_div = a_res[16]
        for j in range(i+1,len(cluster)):
            b_res = cluster[j]
            b_start = b_res[5]
            b_end = b_res[6]
            b_score = b_res[0]
            if ( b_start == 0 and b_end == 0 ):
                continue
            b_len = b_end - b_start + 1
            b_div = b_res[16]
            max_start = max( a_start, b_start)
            min_end = min( a_end, b_end )
            overlap = min_end - max_start + 1
            if ( overlap > 0 ):
               if ( a_len == overlap or b_len == overlap ):
                   # Containment (a in b or b in a ), remove lower
                   # score
                   cluster[j][5] = 0
                   cluster[j][6] = 0
               else:
                   # Overlap
                   if ( a_start < b_start ):
                       # Trim left side
                       cluster[j][5] += overlap
                   else:
                       # Trim right side
                       cluster[j][6] -= overlap



def resolve_using_longer_element( cluster ):
    """
    resolve_using_longer_element( cluster )

    Resolve overlaps by prioritizing longer elements
    over shorter ones. [ greedy ]

    e.g. Resolve the following cluster of overlaps
           ----------------------
                    --------
                          ----------
    as:
           ----------------------
                                 ---
    """
    prev_cres = None
    # Sort by length descending
    cluster = sorted(cluster, key=lambda annot: annot[6]-annot[5], reverse=True)
    for i in range(len(cluster)):
        a_res = cluster[i]
        a_start = a_res[5]
        a_end = a_res[6]
        a_len = a_end - a_start + 1
        if ( a_start == 0 and a_end == 0 ):
            continue
        for j in range(i+1,len(cluster)):
            b_res = cluster[j]
            b_start = b_res[5]
            b_end = b_res[6]
            b_score = b_res[0]
            if ( b_start == 0 and b_end == 0 ):
                continue
            b_len = b_end - b_start + 1
            max_start = max( a_start, b_start)
            min_end = min( a_end, b_end )
            overlap = min_end - max_start + 1
            if ( overlap > 0 ):
               if ( a_len == overlap or b_len == overlap ):
                   # Containment (a in b or b in a ), remove lower
                   # score
                   cluster[j][5] = 0
                   cluster[j][6] = 0
               else:
                   # Overlap
                   if ( a_start < b_start ):
                       # Trim left side
                       cluster[j][5] += overlap
                   else:
                       # Trim right side
                       cluster[j][6] -= overlap


def resolve_using_higher_score( cluster ):
    """
    resolve_using_higher_score( cluster )

    Resolve overlaps by prioritizing higher scoring
    alignments.

    e.g. Resolve the following cluster of overlaps
           -----300--------------
                    --400---
                          ----200---
    as:
                    --400---
                            --200---

    """
    prev_cres = None
    # Sort by score descending
    cluster = sorted(cluster, key=itemgetter(0), reverse=True)
    for i in range(len(cluster)):
        a_res = cluster[i]
        a_start = a_res[5]
        a_end = a_res[6]
        a_len = a_end - a_start + 1
        a_score = a_res[0]
        if ( a_start == 0 and a_end == 0 ):
            continue
        for j in range(i+1,len(cluster)):
            b_res = cluster[j]
            b_start = b_res[5]
            b_end = b_res[6]
            b_score = b_res[0]
            if ( b_start == 0 and b_end == 0 ):
                continue
            b_len = b_end - b_start + 1
            max_start = max( a_start, b_start)
            min_end = min( a_end, b_end )
            overlap = min_end - max_start + 1
            if ( overlap > 0 ):
               if ( a_len == overlap or b_len == overlap ):
                   # Containment (a in b or b in a ), remove lower
                   # score
                   cluster[j][5] = 0
                   cluster[j][6] = 0
               else:
                   # Overlap
                   if ( a_start < b_start ):
                       # Trim left side
                       cluster[j][5] += overlap
                   else:
                       # Trim right side
                       cluster[j][6] -= overlap


#
# main subroutine ( protected from import execution )
#
def main(*args):
    #
    # Options processing
    #
    #   There are two ways to document usage/command line
    #   arguments using this boilerplate.  The intended way
    #   is to document using docstrings at the top of the
    #   script.  This way the pydoc docs match the output
    #   produced by '-h' or '--help' using the argparse
    #   custom action class ( _CustomUsageAction ) defined
    #   below.  If you want the paired-down argparse default
    #   instead, simply remove the "add_help=False" argument
    #   to the argparse constructor below and comment out
    #   the add_argment('-h', ...) line.
    #
    class _CustomUsageAction(argparse.Action):
        """
        _CustomUsageAction() - Class to call our _usage function
        """
        def __init__(self, option_strings, dest, default=False, required=False, help=None):
            super(_CustomUsageAction, self).__init__(
                      option_strings=option_strings, dest=dest,
                      nargs=0, const=True, default=default,
                      required=required, help=help)
        def __call__(self, parser, args, values, option_string=None):
            _usage()

    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('rm_file', metavar='<*.out> or <*.align>')
    parser.add_argument('-h', '--help', action=_CustomUsageAction )
    parser.add_argument("-l", "--log_level", default="INFO")
    parser.add_argument('-d', '--out_dir')
    parser.add_argument('-sp', '--split', type=str, help='Split into files based on name, family, class? This is optional.')
    parser.add_argument('-p', '--out_prefix', type=str, help='Prefix to use for output file - default is first field of input filename')
    parser.add_argument('-m', '--min_length', type=int, help='Minimum size of hit to include in sorted file')
    parser.add_argument('-n', '--min_hit_num', type=int, help='Minimum number of hits in file before being created. Only implemented if --split option is invoked. Optional.')
    parser.add_argument('-dmax', '--max_divergence', type=float, help='Maximum divergence allowed in output file.')
    parser.add_argument('-dmin', '--min_divergence', type=float, help='Minimum divergence allowed in output file.')
    parser.add_argument('-s', '--sort_criterion', type=str, help='Sort criterion, i.e. size, name, family, class, size, or divergence (diverge), etc.')
    parser.add_argument("-o", "--ovlp_resolution", type=str, help='Options are higher_score, longer_element, and lower_divergence. Optional')

    args = parser.parse_args()

    # Setup logging and script timing
    logging.basicConfig(format='')
    logging.getLogger().setLevel(getattr(logging, args.log_level.upper()))
    start_time = time.time()

    LOGGER.info("#\n# RM2Bed.py\n#")

    if ( not os.path.exists(args.rm_file) ):
        raise Exception("File " + args.rm_file + " is missing.")
    LOGGER.info("Data File: " + args.rm_file)

    file_prefix = ""
    if ( args.out_prefix ):
        file_prefix = args.out_prefix
    else:
        fpath, fname = os.path.split(args.rm_file)
        fname = re.sub('\.out$', '', fname, flags=re.IGNORECASE)
        fname = re.sub('\.align$', '', fname, flags=re.IGNORECASE)
        file_prefix = fname

    if ( args.out_dir ):
        if ( not os.path.exists(args.out_dir) ):
            raise Exception("Directory " + args.out_dir + " does not exist.")
        LOGGER.info("Output Directory: " + args.out_dir)
        file_prefix = os.path.join(args.out_dir,file_prefix)
    else:
        LOGGER.info("Output Directory: .")

    # Check ovlp_resolution keywords
    if ( args.ovlp_resolution and args.ovlp_resolution != 'higher_score' and
         args.ovlp_resolution != 'longer_element' and
         args.ovlp_resolution != 'lower_divergence' ):
        raise Exception("--ovlp_resolution keyword '" + args.ovlp_resolution + \
                        "' not recognized.  Must be either 'higher_score', " + \
                        "'longer_element', or 'lower_divergence'")

    ##
    ## Read in file and correct identifiers
    ##   RepeatMasker uses the ID field to join fragments related
    ##   to the same intregration event. In some cases ( notably
    ##   in cluster environments ) it is easier to breakup a genome
    ##   into pieces and run them through RepeatMasker individually
    ##   than it is to run them all at once.  The problem with this
    ##   is that each run restarts ID numbering at 1, making the
    ##   combined *.out, *.align contain redundant IDs.  The code
    ##   below detects changes in the ID number and corrects them.
    ##
    o_rm_file = openOptGzipFile(args.rm_file, modes='r')
    summary_line_RE = re.compile('^\s*\d+\s+\d+\.\d+\s+\d+\.\d+')
    results = []
    cmax_id = 0
    last_rm_id = 0
    new_ids = {}
    last_query_seq = None
    flds = []
    concat_results_detected = 0
    for line in o_rm_file:
        line = line.rstrip()
        # Is this a *.align or *.out alignment summary line?
        if ( summary_line_RE.match(line) ):
            if ( flds ):
                results.append(flds)
            flds = line.split()
            # Out File  :  Always 15 or 16 fields.  The 8th fields is either "C" or "+"
            # Align File:  Forward strand has an empty 8th field while reverse strand
            #              hits have "C" in the eighth field.  The "*" overlap
            #              flag adds the 16th field when it exists.
            if ( len(flds) == 15 or ( len(flds) == 16 and flds[15] == '*' )):
              if ( len(flds) == 16 ):
                del flds[15]
              if ( flds[8] == 'C' ):
                  flds[8] = '-'
              elif ( flds[8] == '+' ):
                  flds[8] = '+'
              else:
                raise Exception("Orientation of RepeatMasker line is unexpected: " + flds[8] )
            elif ( len(flds) == 14 ):
              flds.insert(8,'+')
            else:
              raise Exception("Field count of RepeatMasker line is unexpected: " + str(len(flds)) )
            #
            # Alignment files do not breakup RM identifiers name#type/class
            # into two fields like the *.out files do.  Here we throw out the
            # *.align RM stage identifer column ('m_b#s#i#' or 'c_b#s#i#' )
            # and replace it with the type/class broken out from the combined
            # id.  In this way the datastructure should be the same for both
            # *.out and *.align files.
            if ( 'm_b' in flds[13] or 'c_b' in flds[13] ):
                del flds[13]
                # Alignment file
                if ( '#' in flds[9] ):
                    # name#class/subclass form
                    name, classification = flds[9].split("#")
                    flds[9] = name
                    flds.insert(10,classification)
                else:
                    # class/subclass are not defined
                    flds.insert(10,'unknown')

            # Now breakup the class/subclass into their own columns
            if ( '/' in flds[10] ):
                rmclass, rmsubclass = flds[10].split("/")
                flds[10] = rmclass
                flds.insert(11,rmsubclass)
            else:
                flds.insert(11,'unknown')

            # Fix ID numbers
            #  Two ways this can renumber IDs.  First if the sequence changes
            #  it cannot join fragments between sequences so this can be used
            #  as a natural ID boundary.  The second way we keep the IDs unique
            #  is to detect a fall of over 50 in the ID value coinciding with
            #  a startover of the ID magnitude.
            query_seq = flds[4]
            if ( len(flds) < 16 ):
                print (line)
            rm_id = int(flds[15])
            if ( query_seq != last_query_seq or
                ( rm_id < last_rm_id - 50 and rm_id < 3 ) ):
                if ( query_seq == last_query_seq ):
                    concat_results_detected = 1
                new_ids = {}
            last_rm_id = rm_id
            last_query_seq = query_seq
            if ( rm_id in new_ids ):
                flds[15] = new_ids[rm_id]
            else:
               cmax_id += 1
               new_ids[rm_id] = cmax_id
               flds[15] = cmax_id

            # Convert integer/float fields to native types
            #   -- Do not convert family coordinates as they are not used
            flds[0] = int(flds[0])   # score
            flds[1] = float(flds[1]) # pct mismatch
            flds[2] = float(flds[2]) # pct deletions
            flds[3] = float(flds[3]) # pct insertions
            flds[5] = int(flds[5])   # query start
            flds[6] = int(flds[6])   # query end
            # query remaininig
            flds[7] = int(flds[7].replace('(','').replace(')',''))

            # Finally, create a default divergence column col16
            flds.append(-1.0);

        elif( line.startswith("Kimura") ):
            kDiv = float(line.split('= ')[1])
            flds[16] = kDiv
    if ( flds ):
        results.append(flds)
    o_rm_file.close()
    results = sorted(results, key=itemgetter(5))
    LOGGER.info("Data File Stats:")
    LOGGER.info("   Annotation Lines: " + str(len(results)))
    LOGGER.info("   Insertions (joined frags): " + str(cmax_id))
    if ( concat_results_detected ):
      LOGGER.info("   Info: Concatenated result file detected")

    # Overlap Resolution
    #  - Idea here is that it's faster to use the sorted (by query start)
    #    list produced above to identify clusters of potentially overlapping
    #    annotations first.  Then take that cluster and send it off to a
    #    one of several routines that implement various resolution rules.
    if ( args.ovlp_resolution ):
        LOGGER.info("Overlap Resolution:")
        LOGGER.info("   Method: " + args.ovlp_resolution)
        last_query_seq = None
        max_query_end = 0
        cluster = []
        result_idx = -1
        for result in results:
            result_idx += 1
            query_seq = result[4]
            query_start = result[5]
            query_end = result[6]
            # Overlap detection
            if ( query_seq != last_query_seq or
                 query_start > max_query_end ):
                if ( len(cluster) > 1 ):
                    if ( args.ovlp_resolution == 'higher_score' ):
                        resolve_using_higher_score( cluster )
                    elif ( args.ovlp_resolution == 'longer_element' ):
                        resolve_using_longer_element( cluster )
                    elif ( args.ovlp_resolution == 'lower_divergence' ):
                        resolve_using_lower_divergence( cluster )
                    else:
                        raise Exception("Unknown overlap resolution keyword: " \
                                        + args.ovlp_resolution )
                max_query_end = 0
                cluster = []
            cluster.append(result)
            if ( query_end > max_query_end ):
                max_query_end = query_end
            last_query_seq = query_seq
        # Trailing case
        if ( cluster and len(cluster) > 1 ):
               if ( args.ovlp_resolution == 'higher_score' ):
                   resolve_using_higher_score( cluster )
               elif ( args.ovlp_resolution == 'longer_element' ):
                   resolve_using_longer_element( cluster )
               elif ( args.ovlp_resolution == 'lower_divergence' ):
                   resolve_using_lower_divergence( cluster )
               elif ( args.ovlp_resolution == 'first_element' ):
                   resolve_using_first_element( cluster )
               elif ( args.ovlp_resolution == 'split_between' ):
                   resolve_using_split_between( cluster )
               else:
                   raise Exception("Unknown overlap resolution keyword: " \
                                   + args.ovlp_resolution )
    else:
        LOGGER.info("Overlap Resolution: Keep overlapping annotations")

    # Create a pandas dataframe ... probably not necessary, but it does
    # provide some useful filtering functionality
    annot_dataframe = pd.DataFrame(results, columns=['score',
        'pct_mismatches', 'pct_deletions', 'pct_insertions', 'chrom', 'start',
        'stop', 'rem', 'strand', 'family', 'class', 'subclass', 'unused1',
        'unused2', 'unused3', 'linkage_id', 'diverge'])

    if ( args.ovlp_resolution ):
        # Filter out deleted overlapping annotations.  They are currently
        # marked with query_start = 0 and query_end = 0
        annot_dataframe = annot_dataframe[(annot_dataframe['start'] != 0) & \
                                          (annot_dataframe['stop'] != 0)]
        LOGGER.info("   Remaining annotations: " + str(len(annot_dataframe.index)))

    # BED File coordinates are zero-based, half-open.  RepeatMasker
    # output is 1-based, fully-closed.  Convert to BED coordinates
    annot_dataframe['start'] = annot_dataframe['start'].subtract(1)

    # Size calculated from BED style coordinates
    annot_dataframe['size'] = annot_dataframe['stop'].subtract(annot_dataframe['start'])

    # Columns used for BED output
    #  Field         Desc
    #  ----------    -----------------------------------
    #  sequence      Input sequence
    #  start         Start position 0-based
    #  end           End position 0-based, half-open
    #  family        TE Family Name
    #  score         Raw score from RepeatMasker
    #  orientation   "+"/"-" for forward/reverse strand
    #  class         RepeatMasker Class
    #  subclass      RepeatMasker subclass or "undefined"
    #  divergence    Kimura divergence from *.align file
    #  linkage_id    ID column from RepeatMasker output
    # Simple repeats do not have a divergence calculated
    # for them.  Currently this is marked with the sentinel
    # '-1.0'.
    annot_dataframe = annot_dataframe[['chrom', 'start', 'stop', 'family', 'size', \
                           'strand', 'class', 'subclass', 'diverge', \
                           'linkage_id']]

    # Filter by minimum length
    if ( args.min_length ):
        annot_dataframe = annot_dataframe[annot_dataframe['size'] >= args.min_length]

    # Sort main output if asked.
    if ( args.sort_criterion ):
        LOGGER.info("Sorting By:" + args.sort_criterion)
        if args.sort_criterion in ['family', 'class', 'subclass']:
            annot_dataframe = annot_dataframe.sort_values([args.sort_criterion])
        elif args.sort_criterion in ['size']:
            annot_dataframe = annot_dataframe.sort_values([args.sort_criterion], ascending = [0])
        elif args.sort_criterion in ['diverge']:
            annot_dataframe = annot_dataframe.sort_values([args.sort_criterion], ascending = [1])
        else:
            raise Exception("Invalid sort criterion: " + args.sort_criterion + \
                  ".  Choices are size, family, class, subclass, or " + \
                  "diverge.")
        #log

    # Announce max hits limit
    #if HITS is not None:
    #    print('Will only output files with at least ' + str(HITS) + ' hits.')

    if ( args.min_length or args.max_divergence
         or args.min_divergence ):
        LOGGER.info("Filtering By:")

    # Apply min length filter if requested
    if ( args.min_length ):
        before_cnt = len(annot_dataframe.index)
        annot_dataframe = annot_dataframe[annot_dataframe['size'] >= args.min_length]
        cnt_removed = before_cnt - len(annot_dataframe.index)
        LOGGER.info("   Min Length " + str(args.min_length) + \
                    ": Removed " + str(cnt_removed) + " annotations")

    # Apply max divergence if requested
    if ( args.max_divergence):
        before_cnt = len(annot_dataframe.index)
        annot_dataframe = annot_dataframe[annot_dataframe['diverge'] <= float(args.max_divergence)]
        cnt_removed = before_cnt - len(annot_dataframe.index)
        LOGGER.info("   Max Divergence " + str(args.max_divergence) + \
                    ": Removed " + str(cnt_removed) + " annotations")

    # Apply min divergence if requested
    if ( args.min_divergence ):
        before_cnt = len(annot_dataframe.index)
        annot_dataframe = annot_dataframe[annot_dataframe['diverge'] >= args.min_divergence]
        cnt_removed = before_cnt - len(annot_dataframe.index)
        LOGGER.info("   Min Divergence " + str(args.min_divergence) + \
                    ": Removed " + str(cnt_removed) + " annotations")

    if ( args.min_length or args.max_divergence
         or args.min_divergence ):
        LOGGER.info("   Remaining Annotations: " + str(len(annot_dataframe.index)))

    # Split into files if asked. Also check to see if there is a minumum
    # hit number and act accordingly.
    if ( args.split ):
        LOGGER.info("Split files by: " + args.split)
        if ( args.split in ['name', 'family', 'class'] ):
            clustered = annot_dataframe.sort_values([args.split])
            for split_value in clustered[args.split].unique():
                if args.min_hit_num is not None:
                    clustered_W = clustered[clustered[args.split]==split_value]
                    count_row = clustered_W.shape[0]
                    if count_row >= args.min_hit_num:
                        LOGGER.info("  Creating: " + file_prefix + '_' + \
                                    split_value + '_rm.bed' )
                        clustered_W.to_csv(file_prefix + '_' + split_value + '_rm.bed',
                                           sep='\t', header=False, index=False)
                else:
                    clustered_W = clustered[clustered[args.split]==split_value]
                    LOGGER.info("  Creating: " + file_prefix + '_' + \
                                split_value + '_rm.bed' )
                    clustered_W.to_csv(file_prefix + '_' + split_value + '_rm.bed',
                                       sep='\t', header=False, index=False)
        else:
            print('Splitting options are by name, family, and class.')

    # Write as monolithic file
    LOGGER.info("Creating: " + file_prefix + '_rm.bed' )
    annot_dataframe.to_csv(file_prefix + '_rm.bed', sep='\t', header=False, index=False)


    #
    # Remaining main() code
    #

    end_time = time.time()
    LOGGER.info("Run time: " + str(datetime.timedelta(seconds=end_time-start_time)))


#
# Wrap script functionality in main() to avoid automatic execution
# when imported ( e.g. when help is called on file )
#
if __name__ == '__main__':
    main(*sys.argv)

