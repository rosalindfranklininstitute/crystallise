#
# Copyright (C) 2021 RFI
#
# Author: James Parkhurst
#
# This code is distributed under the GPLv3 license, a copy of
# which is included in the root directory of this package.
#
import json
import mendeleev
import os.path
import urllib.request
import crystallise.data
from html.parser import HTMLParser


class WebElementsHTMLParser(HTMLParser):
    """
    Parse the HTML elements

    """

    def __init__(self):
        super().__init__()
        self._state = None
        self.space_group = None
        self.unit_cell_a = None
        self.unit_cell_b = None
        self.unit_cell_c = None
        self.unit_cell_alpha = None
        self.unit_cell_beta = None
        self.unit_cell_gamma = None

    def handle_data(self, data):
        def parse_default():
            if data.startswith("Space group number:"):
                self._state = "parse_space_group"
            elif data.startswith("Cell parameters:"):
                self._state = "parse_unit_cell"

        def parse_space_group():
            try:
                self.space_group = int(data)
            except:
                pass
            self._state = None

        def parse_unit_cell():
            if data.startswith("a"):
                self._state = "parse_unit_cell_a"
            elif data.startswith("b"):
                self._state = "parse_unit_cell_b"
            elif data.startswith("c"):
                self._state = "parse_unit_cell_c"
            elif data.startswith("α"):
                parse_unit_cell_alpha()
            elif data.startswith("β"):
                parse_unit_cell_beta()
            elif data.startswith("γ"):
                parse_unit_cell_gamma()

        def parse_unit_cell_a():
            tokens = data.split()
            try:
                self.unit_cell_a = float(tokens[1]) / 100
            except:
                pass
            self._state = "parse_unit_cell"

        def parse_unit_cell_b():
            tokens = data.split()
            try:
                self.unit_cell_b = float(tokens[1]) / 100
            except:
                pass
            self._state = "parse_unit_cell"

        def parse_unit_cell_c():
            tokens = data.split()
            try:
                self.unit_cell_c = float(tokens[1]) / 100
            except:
                pass
            self._state = "parse_unit_cell"

        def parse_unit_cell_alpha():
            tokens = data.split()
            try:
                self.unit_cell_alpha = float(tokens[1][:-1])
            except:
                pass
            self._state = "parse_unit_cell"

        def parse_unit_cell_beta():
            tokens = data.split()
            try:
                self.unit_cell_beta = float(tokens[1][:-1])
            except:
                pass
            self._state = "parse_unit_cell"

        def parse_unit_cell_gamma():
            tokens = data.split()
            try:
                self.unit_cell_gamma = float(tokens[1][:-1])
            except:
                pass
            self._state = None

        parse = {
            None: parse_default,
            "parse_space_group": parse_space_group,
            "parse_unit_cell": parse_unit_cell,
            "parse_unit_cell_a": parse_unit_cell_a,
            "parse_unit_cell_b": parse_unit_cell_b,
            "parse_unit_cell_c": parse_unit_cell_c,
            "parse_unit_cell_alpha": parse_unit_cell_alpha,
            "parse_unit_cell_beta": parse_unit_cell_beta,
            "parse_unit_cell_gamma": parse_unit_cell_gamma,
        }[self._state]()


def get_webelements_table():
    """
    Get a table from webelements

    """

    def get_element_list():
        return [e.name.lower() for e in mendeleev.element(list(range(1, 119)))]

    def get_webpage(element):
        url = "https://www.webelements.com/%s/crystal_structure.html" % element
        u2 = urllib.request.urlopen(url)
        return str(u2.read())

    def parse_webpage(page):
        parser = WebElementsHTMLParser()
        parser.feed(page)
        space_group = parser.space_group
        unit_cell = (
            parser.unit_cell_a,
            parser.unit_cell_b,
            parser.unit_cell_c,
            parser.unit_cell_alpha,
            parser.unit_cell_beta,
            parser.unit_cell_gamma,
        )
        if None in unit_cell:
            unit_cell = None
        return space_group, unit_cell

    def save_table(table):
        directory = os.path.dirname(crystallise.data.__file__)
        with open(os.path.join(directory, "elements.json"), "w") as outfile:
            json.dump(table, outfile, indent=2, ensure_ascii=True)

    table = {}
    for Z, element in enumerate(get_element_list(), start=1):
        page = get_webpage(element)
        space_group, unit_cell = parse_webpage(page)
        table[Z] = {
            "space_group": space_group,
            "unit_cell": unit_cell,
        }

    save_table(table)


if __name__ == "__main__":
    get_webelements_table()
