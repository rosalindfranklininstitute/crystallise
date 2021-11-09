#
# Copyright (C) 2021 RFI
#
# Author: James Parkhurst
#
# This code is distributed under the GPLv3 license, a copy of
# which is included in the root directory of this package.
#
import json
import os.path


def element_table():
    """
    Get the element table

    Returns:
        A table of elements where each row contains space group and unit cell

    """
    if not hasattr(element_table, "table"):
        directory = os.path.dirname(__file__)
        with open(os.path.join(directory, "elements.json")) as infile:
            table = json.load(infile)
            table = dict((int(k), v) for k, v in table.items())
            element_table.table = table
    return element_table.table
