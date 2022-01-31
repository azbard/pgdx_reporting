import re
import os


def get_user_input():
    """
    Asks user to either go ahead with latest batch or enter
    batch number. Returns empty string for latest batch or
    "YY-N" (str) for another batch.
    """
    pattern = re.compile(
        r"""       # pattern for entering a batch ID
                \AB?                # "B" or noyhing at start then 
                (?P<year>\d{2})     # a two digit number = year
                -                   # a dash
                (?P<batch>\d+)      # some numbers = batch_num
                \Z                  # the end
                """,
        flags=re.VERBOSE | re.IGNORECASE,
    )
    # the VERBOSE setting ignores whitespace
    # and allows inline comments

    while True:
        user_input = input(
            """
*******************************************

Press Enter to process the latest batch.

Enter a Batch ID (B21-1 or 20-55).

********************************************
"""
        )
        if user_input in ["", "full", "Full"]:
            break
        else:
            m = pattern.match(user_input)
            if not bool(m):
                print("The batch format seems wrong.")
                continue
            else:
                user_input = m.groupdict()
                break

    return user_input


def latest_batch_directory(parent_dir=None):
    """
    Takes a directory (defaults to current) and returns the folder
    name for the latest batch with format BYY-N(NN...)
    Returns "None" (str) if there are no batch folders.
    """

    pattern = re.compile(
        r"""
            \A B                # "B" at start then 
            (?P<year>\d{2})     # a two digit number = year
            -                   # a dash
            (?P<batch>\d+)      # some numbers = batch_num
            \D                  # a non-number
            .+\Z                # any other stuff at the end
            """,
        re.VERBOSE,
    )  # the VERBOSE setting ignores whitespace
    # and allows inline comments

    # initialize
    highest_year = int(0)
    highest_batch = int(0)
    candidate = "None"

    # go through directory
    for dir in os.listdir(parent_dir):

        # check if it's a folder
        if os.path.isdir(os.path.join(parent_dir, dir)):

            # look for the pattern
            m = pattern.search(dir)

            # if you find it
            if bool(m):

                # extract year and batch numbers
                year = int(m.group("year"))
                batch = int(m.group("batch"))

                # if year is higher than highest so far
                # then this is our current candidate
                # for latest batch
                if year > highest_year:
                    highest_year = year
                    highest_batch = batch
                    candidate = dir

                # otherwise just check if batch number
                # is higher
                elif (year == highest_year) and (batch > highest_batch):
                    highest_batch = batch
                    candidate = dir

    return candidate


def get_dir_from_batch(par_dir, batch_dict):
    """
    Takes a directory and a dict{'year':str, 'batch':str} and
    returns the folder name for the associated batch.
    Returns "None" (str) if match not found.
    """

    pattern = re.compile(
        r"""
            \A B                # "B" at start then 
            (?P<year>\d{2})     # a two digit number = year
            -                   # a dash
            (?P<batch>\d+)      # some numbers = batch_num
            \D                  # a non-number
            .+\Z                # any other stuff at the end
            """,
        re.VERBOSE,
    )  # the VERBOSE setting ignores whitespace
    # and allows inline comments

    candidate = "None"

    # go through directory
    for dir in os.listdir(par_dir):

        # check if it's a folder
        if os.path.isdir(os.path.join(par_dir, dir)):

            # look for the pattern
            m = pattern.search(dir)

            # if you find it
            if bool(m):

                # extract year and batch numbers
                year = m.group("year")
                batch = m.group("batch")

                if year == batch_dict["year"] and batch == batch_dict["batch"]:
                    candidate = dir
                    break

    return candidate


def check_batch_already_run(batch_dir):
    """
    Checks if processing already complete.
    """

    igv_dir = os.path.join(batch_dir, "IGV")
    pathologist_dir = os.path.join(batch_dir, "All_Pathologist")
    wiki_dir = os.path.join(batch_dir, "Wiki")

    new_dirs = [wiki_dir, pathologist_dir, igv_dir]

    for dir in new_dirs:
        if os.path.exists(dir):
            return True

    return False


def check_batch_name(to_check: str):
    """
    Asks user to either go ahead with latest batch or enter
    batch number. Returns empty string for latest batch or
    "YY-N" (str) for another batch.
    """
    pattern = re.compile(
        r"""       # pattern for entering a batch ID
                \AB?                # "B" or noyhing at start then 
                (?P<year>\d{2})     # a two digit number = year
                -                   # a dash
                (?P<batch>\d+)      # some numbers = batch_num
                \Z                  # the end
                """,
        flags=re.VERBOSE | re.IGNORECASE,
    )
    # the VERBOSE setting ignores whitespace
    # and allows inline comments

    m = pattern.match(to_check)

    if bool(m):
        return m.groupdict()
    else:
        return None
