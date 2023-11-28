import re
import h5py
from famdb_globals import (
    SOUNDEX_LOOKUP,
    GROUP_FAMILIES,
    dfam_acc_pat,
)
from famdb_helper_classes import Family


def accession_bin(acc):
    """Maps an accession (Dfam or otherwise) into apropriate bins (groups) in HDF5"""
    dfam_match = dfam_acc_pat.match(acc)
    if dfam_match:
        path = (
            GROUP_FAMILIES
            + "/"
            + dfam_match.group(1)
            + "/"
            + dfam_match.group(2)
            + "/"
            + dfam_match.group(3)
            + "/"
            + dfam_match.group(4)
        )
    else:
        path = GROUP_FAMILIES + "/Aux/" + acc[0:2].lower()
    return path


def get_family(entry):
    if not entry:
        return None

    family = Family()

    # Read the family attributes and data
    for k in entry.attrs:
        value = entry.attrs[k]
        setattr(family, k, value)

    return family


def families_iterator(g, prefix=""):
    for key, item in g.items():
        path = "{}/{}".format(prefix, key)
        if isinstance(item, h5py.Dataset):  # test for dataset
            yield (key)
        elif isinstance(item, h5py.Group):  # test for group (go down)
            yield from families_iterator(item, path)


# Filter methods --------------------------------------------------------------------------
def filter_name(family, name):
    """Returns True if the family's name begins with 'name'."""

    if family.attrs.get("name"):
        if family.attrs["name"].lower().startswith(name):
            return True

    return False


def filter_search_stages(family, stages):
    """Returns True if the family belongs to a search stage in 'stages'."""
    if family.attrs.get("search_stages"):
        sstages = (ss.strip() for ss in family.attrs["search_stages"].split(","))
        for family_ss in sstages:
            if family_ss in stages:
                return True

    return False


def filter_repeat_type(family, rtype):
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


def filter_curated(accession, curated):
    """
    Returns True if the family's curatedness is the same as 'curated'. In
    other words, 'curated=True' includes only curated familes and
    'curated=False' includes only uncurated families.

    Families are currently assumed to be curated unless their name is of the
    form DR<9 digits>.

    TODO: perhaps this should be a dedicated 'curated' boolean field on Family
    """

    is_curated = not (
        accession.startswith("DR")
        and len(accession) == 11
        and all((c >= "0" and c <= "9" for c in accession[2:]))
    )

    return is_curated == curated


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
        prev = codes[i - 1]

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
        coding += "0"

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
