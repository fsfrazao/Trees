import pytest
import trees_ibm
import os
import json
import tables

@pytest.fixture
def seedbank():
    seedbank={'FT1':[(2.5,2.5),(3.0,3.0)]}
    return seedbank

@pytest.fixture
def seedbank_on_file(seedbank):
    filename='test_seedbank_on_file.json'
    with open(filename, 'w') as f:
        json.dump(seedbank, f)
    yield filename
    os.remove(filename)


@pytest.fixture
def model_status_on_file():
    model_status = {'FT1': [{"id": 0, "position": [3.5, 3.5], "DBH": 8, "age": 0},
         {"id": 1, "position": [0.5, 1.0], "DBH": 23, "age": 10}],
        'FT2':[{"id": 2, "position": [2.5, 2.5], "DBH": 40, "age": 25},
         {"id": 3, "position": [1.5, 2.0], "DBH": 15, "age": 5}]}
    filename = "test_model_status.json"
    with open(filename, 'w') as f:
        json.dump(model_status, f)
    yield filename
    os.remove(filename)

@pytest.fixture
def database():
    database_name = "test_create_database.h5"
    simple_world.sim_number=1
    simple_world.create_HDF_database(database_name)
    yield database_name
    os.remove(database_name)


def test_load_seedbank(simple_world, seedbank, seedbank_on_file):
    simple_world.seedbank_from_file(seedbank_on_file)
    assert simple_world.topology.seedbank == seedbank


def test_save_seed_bank(simple_world, seedbank):
    filename = "test_save_seed_bank.json"
    simple_world.seedbank_to_file(filename)
    assert os.path.isfile(filename)
    os.remove(filename)


def test_load_model_status(world_with_2_pfts, model_status_on_file):
    world_with_2_pfts.model_status_from_file(model_status_on_file)

    assert trees_ibm.Tree.Instances[0].Ftype == "FT1"
    assert trees_ibm.Tree.Instances[1].Ftype == "FT1"
    assert trees_ibm.Tree.Instances[2].Ftype == "FT2"
    assert trees_ibm.Tree.Instances[3].Ftype == "FT2"

    assert trees_ibm.Tree.Instances[0].age == 0
    assert trees_ibm.Tree.Instances[1].age == 10
    assert trees_ibm.Tree.Instances[2].age == 25
    assert trees_ibm.Tree.Instances[3].age == 5

    assert trees_ibm.Tree.Instances[0].position == (3.5, 3.5)
    assert trees_ibm.Tree.Instances[1].position == (0.5, 1.0)
    assert trees_ibm.Tree.Instances[2].position == (2.5, 2.5)
    assert trees_ibm.Tree.Instances[3].position == (1.5, 2.0)

    assert trees_ibm.Tree.Instances[0].DBH == 8
    assert trees_ibm.Tree.Instances[1].DBH == 23
    assert trees_ibm.Tree.Instances[2].DBH == 40
    assert trees_ibm.Tree.Instances[3].DBH == 15



def test_save_model_status(world_with_2_pfts, model_status_on_file):
    world_with_2_pfts.model_status_from_file(model_status_on_file)
    filename = "test_save_model_status.json"
    world_with_2_pfts.model_status_to_file(filename)

    assert os.path.isfile(filename)
    with open(model_status_on_file) as original_model_status:
        with open(filename) as new_model_status:
            assert original_model_status.read() == new_model_status.read()
    os.remove(filename)


def test_create_database(simple_world):
    database_name = "test_create_database.h5"
    simple_world.sim_number=1
    simple_world.create_HDF_database(database_name)

    assert os.path.isfile(database_name)
    h5 = tables.open_file(database_name)
    h5.get_node("/sim_1")
    h5.get_node("/sim_1/trees")
    h5.get_node("/sim_1/trees/ind_lvl")
    h5.get_node("/sim_1/trees/ind_lvl/Ind")
    h5.get_node("/sim_1/trees/sys_lvl")
    h5.get_node("/sim_1/trees/sys_lvl/NEE")
    h5.get_node("/sim_1/trees/sys_lvl/Pop")
    h5.get_node("/sim_1/trees/sys_lvl/Stocks")
    os.remove(database_name)

def test_write_NEE_data_to_db(world_with_1_tree):
    world_with_1_tree.tsdead_A = 0.1
    world_with_1_tree.tsslow_A = 0.2
    world_with_1_tree.tsfast_A = 0.3
    world_with_1_tree.tsdead_Sslow = 0.4
    world_with_1_tree.tsdead_Sfast = 0.5
    world_with_1_tree.Sdead["current"] = 1.0
    world_with_1_tree.Sdead["previous"] = 1.2
    world_with_1_tree.Sslow["current"] = 1.5
    world_with_1_tree.Sfast["current"] = 2.0
    world_with_1_tree.Cgpp = 2.5
    world_with_1_tree.Cr = 3.0
    world_with_1_tree.Smort = 3.5
    nee = 4.0

    world_with_1_tree.NEE_data_to_db(nee)

    NEE_table = world_with_1_tree.db.get_node(
        "/sim_1/trees/sys_lvl/NEE")
    nee_r = NEE_table[0]

    assert nee_r['dead_wood'] == pytest.approx(0.1)
    assert nee_r['soil_slow'] == pytest.approx(0.3)
    assert nee_r['soil_fast'] == pytest.approx(0.6)
    assert nee_r['gpp'] == pytest.approx(2.5)
    assert nee_r['living_trees'] == pytest.approx(3.0)
    assert nee_r['smort'] == pytest.approx(3.5)
    assert nee_r['t_deadwood_sslow'] == pytest.approx(0.48)
    assert nee_r['t_sslow_sfast'] == pytest.approx(0.6)
    assert nee_r['total_emissions'] == pytest.approx(-1.5)
    assert nee_r['nee'] == pytest.approx(4.0)





def test_write_Stocks_data_to_db(world_with_2_populations):
    AGB = trees_ibm.Tree.Total_AGB()
    world_with_2_populations.Total_AGB = AGB * 0.44
    world_with_2_populations.Sdead["current"] = 1.2
    world_with_2_populations.Sslow["current"] = 1.5
    world_with_2_populations.Sfast["current"] = 2.0

    world_with_2_populations.stocks_data_to_db()

    Stocks_table = world_with_2_populations.db.get_node(
        "/sim_1/trees/sys_lvl/Stocks")
    stocks_r = Stocks_table[0]

    assert stocks_r['agb'] == pytest.approx(AGB * 0.44)
    assert stocks_r['dead_wood'] == pytest.approx(1.2)
    assert stocks_r['soil_slow'] == pytest.approx(1.5)
    assert stocks_r['soil_fast'] == pytest.approx(2.0)

def test_write_Pop_data_to_db(world_with_2_populations):
    world_with_2_populations.pop_data_to_db()

    Pop_table = world_with_2_populations.db.get_node(
        "/sim_1/trees/sys_lvl/Pop")
    pop_r = Pop_table[0]

    assert pop_r["FT1"] == 10
    assert pop_r["FT1"] == 10

def test_write_Ind_data_to_db(world_with_1_tree):
    world_with_1_tree.ind_data_to_db(1)

    Ind_table = world_with_1_tree.db.get_node(
        "/sim_1/trees/ind_lvl/Ind")
    ind_r = Ind_table[0]

    assert ind_r["step"] == 1
    assert ind_r["ind_id"] == 0
    assert ind_r["pos_x"] == 2.5
    assert ind_r["pos_y"] == 2.0
    assert ind_r["age"] == 0
    assert ind_r["dbh"] == 10
    assert ind_r["agb"] == pytest.approx(0.041243851)
    assert ind_r["height"] == pytest.approx(13.137537)
    assert ind_r["crown_area"] == pytest.approx(6.4772758)
    assert ind_r["FT"] == b"FT1"


def test_create_agents_default(simple_world,pft1):
    FT1 = trees_ibm.trees_ibm.Tree.TreeFactory(new_cls_name="FT1",
    new_parameters=pft1)
    simple_world.create_agents(FT=FT1, n=3, dbh=2, world=simple_world)
    
    assert len(trees_ibm.Tree.Instances.items()) == 3
    
    trees_ibm.Tree.Instances.clear()
    trees_ibm.Tree.ID=0

def test_create_agents_defined_positions(simple_world,pft1):
    FT1 = trees_ibm.trees_ibm.Tree.TreeFactory(new_cls_name="FT1",
    new_parameters=pft1)
    in_positions=[(1.0,1.5),(2.0,2.5)]
    simple_world.create_agents(FT=FT1, n=2, pos=in_positions,  dbh=2, world=simple_world)
    
    assert len(trees_ibm.Tree.Instances.items()) == 2
    out_positions=[t.position for t in trees_ibm.Tree.Instances.values()]
    assert out_positions == in_positions
    
    trees_ibm.Tree.Instances.clear()
    trees_ibm.Tree.ID=0

def test_create_agents_defined_ids_ages_dbhs(simple_world,pft1):
    FT1 = trees_ibm.trees_ibm.Tree.TreeFactory(new_cls_name="FT1",
    new_parameters=pft1)
    in_positions=[(2.0,2.5),(3.0,3.5)]
    in_ids=[75,103]
    in_ages=[77,88]
    in_dbhs=[100,90]
    
    
    simple_world.create_agents(FT=FT1, n=2, pos=in_positions, ids=in_ids, ages=in_ages, DBHs=in_dbhs, world=simple_world)
    
    assert len(trees_ibm.Tree.Instances.items()) == 2
    
    out_positions=[t.position for t in trees_ibm.Tree.Instances.values()]
    out_ids=[t.id for t in trees_ibm.Tree.Instances.values()]
    out_ages=[t.age for t in trees_ibm.Tree.Instances.values()]
    out_dbhs=[t.DBH for t in trees_ibm.Tree.Instances.values()]
    assert out_positions == in_positions
    assert out_ids == in_ids
    assert out_ages == in_ages
    assert out_dbhs == in_dbhs
    
    trees_ibm.Tree.Instances.clear()
    trees_ibm.Tree.ID=0