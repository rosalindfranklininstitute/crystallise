#
# Copyright (C) 2021 RFI
#
# Author: James Parkhurst
#
# This code is distributed under the GPLv3 license, a copy of
# which is included in the root directory of this package.
#
import gemmi
import numpy
import os.path
import crystallise.data
from math import sqrt, ceil

__all__ = [
    "CoordType",
    "is_unit_cell_consistent_with_space_group",
    "generate_fractional_atom_positions",
    "generate_crystal_lattice",
    "filter_atoms",
    "generate_atomic_model",
    "write_coordinates",
]


# The coordinate data type
CoordType = [("atomic_number", "i4"), ("x", "f4"), ("y", "f4"), ("z", "f4")]


def get_unit_cell(atomic_number):
    """
    Get the standard unit cell for a given element.

    From https://www.webelements.com/${name}/crystal_structure.html

    Params:
        atomic_number (int): The atomic number

    Returns:
        tuple: The unit cell parameters

    """
    return crystallise.data.element_table()[atomic_number]["unit_cell"]


def get_space_group(atomic_number):
    """
    Get the standard space group for a given element.

    From https://www.webelements.com/${name}/crystal_struture.html

    Params:
        atomic_number (int): The atomic number

    Returns:
        int: The space group number

    """
    return crystallise.data.element_table()[atomic_number]["space_group"]


def is_unit_cell_consistent_with_space_group(unit_cell, space_group):
    """
    Perform a simple check to see if the crystal system is consistent with the
    given space group.

    Params:
        unit_cell (object): The unit cell object
        space_group (object): The space group object

    Returns:
        bool: Are they consistent

    """

    # Get the crystal system
    crystal_system = space_group.crystal_system()

    # Check floating point values are approximately the same
    def approx(a, b, eps=1e-7):
        return abs(a - b) < eps

    # Check the unit cell against the crystal system
    # fmt: off
    return {
        gemmi.CrystalSystem.Triclinic: lambda cell: True,
        gemmi.CrystalSystem.Monoclinic: lambda cell: (
            approx(cell.alpha, 90)
            and approx(cell.gamma, 90)
        ),
        gemmi.CrystalSystem.Orthorhombic: lambda cell: (
            approx(cell.alpha, 90)
            and approx(cell.beta, 90)
            and approx(cell.gamma, 90)
        ),
        gemmi.CrystalSystem.Tetragonal: lambda cell: (
            approx(cell.a, cell.b)
            and approx(cell.alpha, 90)
            and approx(cell.beta, 90)
            and approx(cell.gamma, 90)
        ),
        gemmi.CrystalSystem.Trigonal: lambda cell: (
            approx(cell.a, cell.b)
            and approx(cell.a, cell.c)
            and approx(cell.alpha, cell.beta)
            and approx(cell.alpha, cell.gamma)
        ),
        gemmi.CrystalSystem.Hexagonal: lambda cell: (
            approx(cell.a, cell.b)
            and approx(cell.alpha, 90)
            and approx(cell.beta, 90)
            and approx(cell.gamma, 120)
        ),
        gemmi.CrystalSystem.Cubic: lambda cell: (
            approx(cell.a, cell.b)
            and approx(cell.a, cell.c)
            and approx(cell.alpha, 90)
            and approx(cell.beta, 90)
            and approx(cell.gamma, 90)
        ),
    }[crystal_system](unit_cell)
    # fmt: on


def generate_fractional_atom_positions(asymmetric_unit, space_group):
    """
    Given the atoms in the ASU and the space group, generate the fractional
    atom positions.

    Params:
        asymmetric_unit (list): The atoms in the asymmetric unit
        space_group (object): The space group object

    Returns:
        array: The fraction atom positions in the unit cell

    """

    # Need to fix this. Currently only handle single atom case
    if asymmetric_unit.shape[0] != 1:
        raise RuntimeError("Currently only handle a single atom")

    # Get the atom positions
    coordinates = []
    for op in space_group.operations():
        for Z, x, y, z in asymmetric_unit:
            x, y, z = op.apply_to_xyz((x, y, z))
            coordinates.append((Z, x, y, z))
    return numpy.array(list(set(coordinates)), dtype=CoordType)


def generate_crystal_lattice(unit_cell, num_unit_cells, fractional_coordinates=None):
    """
    Generate a crystal lattice

    Params:

        unit_cell (object): The unit cell object
        num_unit_cells (tuple): The number of unit cells along each axis
        fractional_coordinates (array): The fractional coordinates within a single unit cell

    Returns:
        list: The list of positions relative to the origin in real space

    """

    # Unpack the number in each dimension
    num_a, num_b, num_c = num_unit_cells

    # Init the fractional coordinates
    if fractional_coordinates is None:
        fractional_coordinates = numpy.array([[0, 0, 0, 0]])

    # Generate the crystal lattice
    shape = num_a * num_b * num_c * fractional_coordinates.shape[0]
    coordinates = numpy.zeros(shape=shape, dtype=CoordType)
    count = 0
    for h in range(num_a):
        for k in range(num_b):
            for l in range(num_c):
                for Z, x, y, z in fractional_coordinates:
                    coordinates[count] = (Z,) + tuple(
                        unit_cell.orthogonalize(gemmi.Fractional(h + x, k + y, l + z))
                    )
                    count += 1

    # Return the coordinates
    return coordinates


def filter_atoms(coordinates, num_atoms=None, morphology="sphere"):
    """
    Filter the atoms so that the crystal has a specific morphology with a given number of atoms

    Params:
        coordinates (array): The atom coordinates
        num_atoms (int): The number of atoms
        morphology (str): The morphology of the crystal

    Returns:
        array: The filtered coordinates

    """

    def filter_atoms_sphere(coordinates, num_atoms):

        # Get the centre of mass
        x = coordinates["x"]
        y = coordinates["y"]
        z = coordinates["z"]
        c = numpy.array([x, y, z]).T
        centre_of_mass = numpy.mean(c, axis=0)

        # Compute all square distances
        sq_distance = numpy.sum((c - centre_of_mass) ** 2, axis=1)

        # Get the selection of the closest n atoms
        index = numpy.argsort(sq_distance)
        return coordinates[index[0:num_atoms]]

    # If the number of atoms is not set then return as is
    if num_atoms is None or morphology is None:
        return coordinates

    # Filter the atoms into the given morphology
    return {"sphere": filter_atoms_sphere}[morphology](coordinates, num_atoms)


def generate_atomic_model(
    unit_cell_parameters,
    space_group_name,
    num_unit_cells,
    asymmetric_unit,
    morphology=None,
    num_atoms=None,
):
    """
    Generate the atomic model of the crystal

    Params:
        unit_cell_parameters (tuple): The unit cell parameters (1, 3 or 6 elements)
        space_group_name (str): The name of the space group
        num_unit_cells (tuple): The number of unit cells along each axis
        asymmetric_unit (object): The contents of the asymmetric unit (str or array)
        morphology (str): The crystal morphology
        num_atoms (int): The number of atoms if morphology is specified

    Returns:
        array: The atomic model of the crystal

    """

    # Check the contents of the asymmetric unit
    if isinstance(asymmetric_unit, str):

        # Get the atomic number and init the array
        Z = gemmi.Element(asymmetric_unit).atomic_number
        asymmetric_unit = numpy.array(
            [(Z, 0, 0, 0)],
            dtype=CoordType,
        )

        # If we have no unit cell parameters or space group try to get them
        # from the known table of unit cells and space groups for elements
        if unit_cell_parameters is None:
            unit_cell_parameters = get_unit_cell(Z)
        if space_group_name is None:
            space_group_name = get_space_group(Z)

    # Get the space group object
    if isinstance(space_group_name, int):
        space_group = gemmi.SpaceGroup(space_group_name)
    else:
        space_group = gemmi.find_spacegroup_by_name(space_group_name)
    if space_group is None:
        raise RuntimeError("Unrecognised space group: %s" % space_group_name)

    # Check the length of the unit cell parameters list
    if len(unit_cell_parameters) == 1:
        unit_cell_parameters = (
            unit_cell_parameters[0],
            unit_cell_parameters[0],
            unit_cell_parameters[0],
            90,
            90,
            90,
        )
    elif len(unit_cell_parameters) == 3:
        unit_cell_parameters = (
            unit_cell_parameters[0],
            unit_cell_parameters[1],
            unit_cell_parameters[2],
            90,
            90,
            90,
        )
    else:
        assert len(unit_cell_parameters) == 6

    # Create the unit cell object
    unit_cell = gemmi.UnitCell(*unit_cell_parameters)

    # Check the unit cell is consistent with the space group
    assert is_unit_cell_consistent_with_space_group(unit_cell, space_group)

    # Generate the fractional atom positions in the unit cell from the asu
    fractional_coordinates = generate_fractional_atom_positions(
        asymmetric_unit, space_group
    )

    # Compute number of unit cells
    if num_unit_cells is None:
        num_unit_cells = int(
            ceil(
                sqrt(3)
                * (float(num_atoms) / len(fractional_coordinates)) ** (1.0 / 3.0)
            )
        )
    if isinstance(num_unit_cells, int):
        num_unit_cells = (num_unit_cells, num_unit_cells, num_unit_cells)
    else:
        assert len(num_unit_cells) == 3

    # Generate the atomic coordinates
    coordinates = generate_crystal_lattice(
        unit_cell, num_unit_cells, fractional_coordinates=fractional_coordinates
    )

    # Filter the coordinates
    return filter_atoms(coordinates, num_atoms, morphology=morphology)


def write_coordinates(coordinates, filename):
    """
    Write the coordinates to file

    Params:
        coordinates (array): The atomic model coordinates
        filename (str): The filename string (either a .pdb, .cif or .csv file)

    """

    # Get a gemmi structure object
    def get_structure(coordinates):
        residue = gemmi.Residue()
        for Z, x, y, z in coordinates:
            atom = gemmi.Atom()
            atom.element = gemmi.Element(Z)
            atom.pos.x = x
            atom.pos.y = y
            atom.pos.z = z
            atom.occ = 1
            atom.charge = 0
            atom.b_iso = 0
            residue.add_atom(atom)
        chain = gemmi.Chain("1")
        chain.add_residue(residue)
        model = gemmi.Model("1")
        model.add_chain(chain)
        structure = gemmi.Structure()
        structure.add_model(model)
        return structure

    # Write a PDB file
    def write_pdb(coordinates, filename):
        get_structure(coordinates).write_minimal_pdb(filename)

    # Write a CIF file
    def write_cif(coordinates, filename):
        get_structure(coordinates).make_mmcif_document().write_file(filename)

    # Write a CSV file
    def write_csv(coordinates, filename):
        with open(filename, "w") as outfile:
            for Z, x, y, z in coordinates:
                outfile.write("%d, %.3f, %.3f, %.3f\n" % (Z, x, y, z))

    # Write the coordinates to file
    {".pdb": write_pdb, ".cif": write_cif, ".csv": write_csv}[
        os.path.splitext(filename)[1].lower()
    ](coordinates, filename)
