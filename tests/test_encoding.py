import numpy as np
from pepdesign.encoding import (AA, one_hot, physchem, combined_features,
                                composition_vector, chain_adjacency,
                                clean_sequence, PHYSCHEM_DIM)


def test_one_hot_shape_and_content():
    m = one_hot("ACD", 10)
    assert m.shape == (10, len(AA))
    assert m[0].sum() == 1.0 and m[1].sum() == 1.0 and m[9].sum() == 0.0
    assert m[0, 0] == 1.0  # A is index 0


def test_one_hot_ignores_noncanonical():
    # noncanonical letters are dropped before encoding: "AXBZC" -> "AC"
    m = one_hot("AXBZC", 5)
    assert m[0].sum() == 1.0 and m[1].sum() == 1.0
    assert m[2].sum() == 0.0


def test_physchem_values():
    m = physchem("KR", 4)
    assert m.shape == (4, PHYSCHEM_DIM)
    assert m[0, 1] == 1.0   # K charge +1
    assert m[1, 1] == 1.0   # R charge +1


def test_combined_features_dim():
    assert combined_features("ACDEFG", 30).shape == (30, 24)


def test_composition_normalized():
    v = composition_vector("AAAAAA")
    assert abs(v.sum() - 1.0) < 1e-9
    assert v.shape == (len(AA) ** 2,)


def test_chain_adjacency_symmetric_normalized():
    a = chain_adjacency(8, window=2)
    assert np.allclose(a, a.T)
    assert a.shape == (8, 8)
    assert (a.diagonal() > 0).all()  # self loops


def test_clean_sequence_uppercases():
    assert clean_sequence("acdxb") == "ACD"
