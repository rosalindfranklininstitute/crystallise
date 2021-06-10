#
# Copyright (C) 2021 RFI
#
# Author: James Parkhurst
#
# This code is distributed under the GPLv3 license, a copy of
# which is included in the root directory of this package.
#
import argparse
import gemmi
import logging
import crystallise.crystallise


def main(args=None):

    # The command line parser
    parser = argparse.ArgumentParser(
        prog="crystallise", description="Generate atomic models of crystals"
    )

    # Add the command line options
    parser.add_argument(
        "-uc",
        "--unit_cell",
        dest="unit_cell",
        type=lambda s: tuple([float(x) for x in s.split(",")]),
        default=(1, 1, 1, 90, 90, 90),
        help="The unit cell parameters",
    )
    parser.add_argument(
        "-sg",
        "--space_group",
        dest="space_group",
        type=str,
        default="P1",
        help="The space group",
    )
    parser.add_argument(
        "-nu",
        "--num_unit_cells",
        dest="num_unit_cells",
        type=lambda s: tuple([int(x) for x in s.split(",")]),
        default=(1, 1, 1),
        help="The number of unit cells along each axis",
    )
    parser.add_argument(
        "-asu",
        "--asymmetric_unit",
        dest="asymmetric_unit",
        type=str,
        default="C",
        help="The contents of the asymmetric unit",
    )
    parser.add_argument(
        "-mo",
        "--morphology",
        dest="morphology",
        type=str,
        choices=["none", "sphere"],
        default=None,
        help="The morphology of the crystal",
    )
    parser.add_argument(
        "-na",
        "--num_atoms",
        dest="num_atoms",
        type=int,
        default=None,
        help="The number of atoms to filter if given a morphology",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=str,
        default="atoms.pdb",
        help="The output file with the atomic model (PDB/CIF/CSV format)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Set verbose output",
    )

    # Parse some argument lists
    args = parser.parse_args(args=args)

    # Set the logger
    if args.verbose is True:
        level = logging.INFO
    else:
        level = logging.WARN
    logging.basicConfig(level=level, format="%(msg)s")

    # Create the atomic model
    coordinates = crystallise.crystallise.generate_atomic_model(
        args.unit_cell,
        args.space_group,
        args.num_unit_cells,
        args.asymmetric_unit,
        args.morphology,
        args.num_atoms,
    )

    # Save the atomic model
    crystallise.crystallise.write_coordinates(coordinates, args.output)
