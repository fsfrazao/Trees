import pytest
import trees_ibm
import os
import tables


@pytest.fixture
def pft1():
    """ Parameters for one plant functional group. """
    pft_par = dict(
        max_agb=200,
        max_dbh=145,
        h0=3.3,
        h1=0.60,
        cl0=0.80,
        cd0=0.60,
        cd1=0.68,
        cd2=0.0,
        rho=0.55,
        sigma=0.70,
        f0=0.77,
        f1=-0.18,
        l0=2.0,
        l1=0.10,
        Iseed=0.03,
        Nseed=20,
        m=0.5,
        rg=0.25,
        Mb=0.015,
        pmax=2,
        alpha=0.36,
        Dmort=10,
        Dfall=0.45,
        pfall=0.3,
        DdeltaDmax=0.33,
        deltaDmax=0.012,
        m_max=0.12,
        Fdisp=1.0,
        Adisp=0.118,
        Nfruit=30)

    return pft_par


@pytest.fixture
def pft2():
    """ Parameters for one plant functional group. """
    pft_par = dict(
        max_agb=200,
        max_dbh=58,
        h0=4.6,
        h1=0.4,
        cl0=0.80,
        cd0=0.60,
        cd1=0.68,
        cd2=0.0,
        rho=0.55,
        sigma=0.70,
        f0=0.77,
        f1=-0.18,
        l0=2.0,
        l1=0.10,
        Iseed=0.01,
        Nseed=15,
        m=0.5,
        rg=0.25,
        Mb=0.03,
        pmax=3.1,
        alpha=0.28,
        Dmort=10,
        Dfall=0.45,
        pfall=0.3,
        DdeltaDmax=0.34,
        deltaDmax=0.012,
        m_max=0.12,
        Fdisp=1.0,
        Adisp=0.118,
        Nfruit=30)

    return pft_par


@pytest.fixture
def simple_world():
    """ A world without any trees or database"""
    topology = trees_ibm.Tree_Grid(x_max=5, y_max=5, delta_h=0.5,
        patch_area=400, I0=860, lday=12, phi_act=365, k=0.7)
    world = trees_ibm.Tree_World(topology=topology)

    yield world
    del world
    trees_ibm.Tree.Instances.clear()
    trees_ibm.Tree.ID=0


@pytest.fixture
def world_with_1_tree(simple_world, pft1):
    FT1 = trees_ibm.trees_ibm.Tree.TreeFactory(new_cls_name="FT1",
    new_parameters=pft1)
    database_name = "Test_world_with_1_tree.h5"
    simple_world.create_HDF_database(database_name)
    simple_world.db = tables.open_file(database_name,"r+")
    simple_world.create_agents(FT=FT1, pos=[(2.5, 2.0)], n=1, dbh=10, world=simple_world)
    simple_world.initial_update()
    yield simple_world
    os.remove(database_name)
    del simple_world
    trees_ibm.Tree.Instances.clear()
    trees_ibm.Tree.ID=0


@pytest.fixture
def one_tree(world_with_1_tree):
    tree = trees_ibm.Tree.Instances.get(0)
    tree.update()
    yield tree


@pytest.fixture
def world_with_2_pfts(simple_world, pft1, pft2):
    """ World instance with 2 pfts (but no trees)"""
    FT1 = trees_ibm.trees_ibm.Tree.TreeFactory(new_cls_name="FT1",
    new_parameters=pft1)
    FT2 = trees_ibm.trees_ibm.Tree.TreeFactory(new_cls_name="FT2",
    new_parameters=pft2)
    yield simple_world
    del simple_world
    trees_ibm.Tree.Instances.clear()
    trees_ibm.Tree.ID=0
    #TODO: clean DERIVED_TREES


@pytest.fixture
def world_with_2_populations(simple_world, pft1, pft2):
    FT1 = trees_ibm.Tree.TreeFactory(new_cls_name="FT1",
    new_parameters=pft1)
    FT2 = trees_ibm.Tree.TreeFactory(new_cls_name="FT2",
    new_parameters=pft2)
    database_name = "Test_world_with_2_populations.h5"
    simple_world.create_HDF_database(database_name)
    simple_world.db = tables.open_file(database_name,"r+")
    # simple_world.topology.update()
    simple_world.create_agents(FT=FT1, n=10, dbh=10, world=simple_world)
    simple_world.create_agents(FT=FT2, n=10, dbh=10, world=simple_world)
    #simple_world.initial_update()
    yield simple_world
    os.remove(database_name)
    del simple_world
    trees_ibm.Tree.Instances.clear()
    trees_ibm.Tree.ID=0
    #TODO: clean DERIVED_TREES
