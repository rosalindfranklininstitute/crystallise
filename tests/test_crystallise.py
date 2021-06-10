import crystallise
import numpy


def test_generate_atomic_model_1():

    unit_cell = (10,)
    space_group = "P1"
    num_unit_cells = 2
    asymmetric_unit = "C"
    morphology = None
    num_atoms = None

    coordinates = crystallise.generate_atomic_model(
        unit_cell,
        space_group,
        num_unit_cells,
        asymmetric_unit,
        morphology,
        num_atoms,
    )

    assert len(coordinates) == (num_unit_cells) ** 3


def test_generate_atomic_model_2():

    unit_cell = (10, 20, 30)
    space_group = "P1"
    num_unit_cells = (3, 4, 5)
    asymmetric_unit = "C"
    morphology = None
    num_atoms = None

    coordinates = crystallise.generate_atomic_model(
        unit_cell,
        space_group,
        num_unit_cells,
        asymmetric_unit,
        morphology,
        num_atoms,
    )

    assert len(coordinates) == numpy.prod(num_unit_cells)


def test_generate_atomic_model_3():

    unit_cell = (4, 4, 4, 90, 90, 90)
    space_group = "Fm-3m"
    num_unit_cells = (3, 4, 5)
    asymmetric_unit = "Au"
    morphology = None
    num_atoms = None

    coordinates = crystallise.generate_atomic_model(
        unit_cell,
        space_group,
        num_unit_cells,
        asymmetric_unit,
        morphology,
        num_atoms,
    )

    assert len(coordinates) == numpy.prod(num_unit_cells) * 4


def test_generate_atomic_model_4():

    unit_cell = (4, 4, 4, 90, 90, 90)
    space_group = "Fm-3m"
    num_unit_cells = (3, 4, 5)
    asymmetric_unit = "Au"
    morphology = "sphere"
    num_atoms = 10

    coordinates = crystallise.generate_atomic_model(
        unit_cell,
        space_group,
        num_unit_cells,
        asymmetric_unit,
        morphology,
        num_atoms,
    )

    assert len(coordinates) == num_atoms


def test_generate_atomic_model_5():

    unit_cell = (4, 4, 4, 90, 90, 90)
    space_group = "Fm-3m"
    num_unit_cells = None
    asymmetric_unit = "Au"
    morphology = "sphere"
    num_atoms = 10

    coordinates = crystallise.generate_atomic_model(
        unit_cell,
        space_group,
        num_unit_cells,
        asymmetric_unit,
        morphology,
        num_atoms,
    )

    assert len(coordinates) == num_atoms


def test_generate_atomic_model_6():

    unit_cell = None
    space_group = None
    num_unit_cells = None
    asymmetric_unit = "Au"
    morphology = "sphere"
    num_atoms = 10

    coordinates = crystallise.generate_atomic_model(
        unit_cell,
        space_group,
        num_unit_cells,
        asymmetric_unit,
        morphology,
        num_atoms,
    )

    assert len(coordinates) == num_atoms


def test_write_coordinates(tmpdir):

    unit_cell = (4, 4, 4, 90, 90, 90)
    space_group = "Fm-3m"
    num_unit_cells = None
    asymmetric_unit = "Au"
    morphology = "sphere"
    num_atoms = 10

    coordinates = crystallise.generate_atomic_model(
        unit_cell,
        space_group,
        num_unit_cells,
        asymmetric_unit,
        morphology,
        num_atoms,
    )

    assert len(coordinates) == num_atoms

    crystallise.write_coordinates(coordinates, str(tmpdir.join("temp.pdb")))
    assert tmpdir.join("temp.pdb").exists()

    crystallise.write_coordinates(coordinates, str(tmpdir.join("temp.cif")))
    assert tmpdir.join("temp.cif").exists()

    crystallise.write_coordinates(coordinates, str(tmpdir.join("temp.csv")))
    assert tmpdir.join("temp.csv").exists()
