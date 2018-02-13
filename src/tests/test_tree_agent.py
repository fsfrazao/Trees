import pytest
import trees_ibm
import os

def test_create_PFT(pft1,simple_world):
    FT1 = trees_ibm.trees_ibm.Tree.TreeFactory(new_cls_name="FT1",
    new_parameters=pft1)


    assert trees_ibm.trees_ibm.Tree.DERIVED_TREES.get("FT1") == FT1
    #class attributes
    assert FT1.Nseed == 20
    assert FT1.Nfruit == 30
    assert FT1.Fdisp == 1.0
    assert FT1.Adisp == 0.118
    assert FT1.Ftype == "FT1"
    assert FT1.Iseed == 0.03
    assert FT1.h0 == 3.3
    assert FT1.h1 == 0.60
    #instance Attributes
    simple_world.create_agents(FT=FT1, pos=[(2.5, 2.0)], n=1, dbh=10, world = simple_world)
    tree = trees_ibm.Tree.Instances.get(0)

    assert tree.rho == 0.55
    assert tree.cl0 == 0.80
    assert tree.cd0 == 0.60
    assert tree.cd1 == 0.68
    assert tree.cd2 == 0.0
    assert tree.sigma == 0.70
    assert tree.f0 == 0.77
    assert tree.f1 == -0.18
    assert tree.l0 == 2.0
    assert tree.l1 == 0.10
    assert tree.k == 0.7
    assert tree.lday == 12
    assert tree.phi_act == 30
    assert tree.m == 0.5
    assert tree.alpha == 0.36
    assert tree.pmax == 2
    assert tree.Dmax == (145 / 100.)
    assert tree.deltaDmax == 0.012
    assert tree.DdeltaDmax == 0.33
    assert tree.Dmort == 10
    assert tree.Dfall == 0.45
    assert tree.pfall == 0.3
    assert tree.rg == 0.25
    assert tree.m_max == 0.12


@pytest.mark.parametrize("Nfruit,DBH,expected", [
    (30, 2, 0),
    (50, 20, 10),
    (0, 2, 0),
    (30.0, 0, 0)

])
def test_fruits_from_DBH(one_tree, Nfruit, DBH, expected):
    one_tree.Nfruit = Nfruit
    one_tree.DBH = DBH
    assert one_tree.fruits_from_DBH() == expected


@pytest.mark.parametrize("DBH, expected", [
    (2, 5.002),
    (20.0, 19.91),
    (0, 0),

])
def test_H_from_DBH(one_tree, DBH, expected):
    one_tree.DBH = DBH
    assert one_tree.H_from_DBH() == pytest.approx(expected, 0.1)


@pytest.mark.parametrize("H, expected", [
    (2, 1.6),
    (20, 16.0),
    (0, 0)
])
def test_CL_from_H(one_tree, H, expected):
    one_tree.H = H
    assert one_tree.CL_from_H() == pytest.approx(expected, 0.1)


@pytest.mark.parametrize("CD, expected", [
    (2, 3.14),
    (20, 314.15),
    (0, 0)
])
def test_CA_from_CD(one_tree, CD, expected):
    one_tree.CD = CD
    assert one_tree.CA_from_CD() == pytest.approx(expected, 0.1)


@pytest.mark.parametrize("H, DBH, expected", [
    (5, 2,  0.0015),
    (19, 20, 0.6),
    (0, 0, 0)
])
def test_calculate_volume(one_tree, H, DBH, expected):
    one_tree.H = H
    one_tree.DBH = DBH
    assert one_tree.calculate_volume() == pytest.approx(expected, 0.1)


@pytest.mark.parametrize("DBH, expected", [
    (0.5, 1.031),
    (2, 16.497),
    (0, 0)
])
def test_AGB_from_DBH(one_tree, DBH, expected):
    assert one_tree.AGB_from_DBH(DBH) == pytest.approx(expected, 0.1)


@pytest.mark.parametrize("DBH, expected", [
    (0.5, 0.8723),
    (2, 0.6796 )
])
def test_calculate_f(one_tree, DBH, expected):
     one_tree.DBH = DBH
     assert one_tree.calculate_f() == pytest.approx(expected, 0.1)


@pytest.mark.parametrize("DBH, expected", [
    (0.5,1.86),
    (2, 2.14),
    (0, 0)
])
def test_LAI_from_DBH(one_tree, DBH, expected):
    one_tree.DBH = DBH
    assert one_tree.LAI_from_DBH() == pytest.approx(expected, 0.1)


@pytest.mark.parametrize("fruit_stock, n, expected", [
    (30, 5, 25),
    (30, 0, 30),
    (0, 1, -1)
])
def test_decrement_fruit_stock(one_tree, fruit_stock, n, expected):
    one_tree.available_fruits = fruit_stock
    one_tree.decrement_fruit_stock(n=n)
    assert one_tree.available_fruits == pytest.approx(expected, 0.1)
