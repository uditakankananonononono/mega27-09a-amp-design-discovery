import numpy as np
import torch
from pepdesign.encoding import combined_features, physchem, chain_adjacency
from pepdesign.models.cnn import PepCNN, count_parameters
from pepdesign.models.gnn import PepGNN, graph_inputs
from pepdesign.data.fixtures import AMP_FIXTURE, NON_AMP_FIXTURE

SEQS = AMP_FIXTURE + NON_AMP_FIXTURE
LABELS = [1] * len(AMP_FIXTURE) + [0] * len(NON_AMP_FIXTURE)
MAX_LEN = 40


def cnn_inputs(seqs, max_len):
    arr = np.stack([combined_features(s, max_len) for s in seqs])
    return (torch.from_numpy(arr),)


def test_cnn_forward_shape():
    model = PepCNN()
    out = model(cnn_inputs(SEQS, MAX_LEN)[0])
    assert out.shape == (len(SEQS),)
    assert count_parameters(model) > 1000


def test_gnn_forward_shape():
    model = PepGNN()
    feats, adj = graph_inputs(SEQS, MAX_LEN)
    assert feats.shape == (len(SEQS), MAX_LEN, 4)
    out = model(feats, adj)
    assert out.shape == (len(SEQS),)


def test_cnn_learns_fixture_signal():
    from pepdesign.models.ensemble import train_binary, predict_scores
    model = PepCNN()
    train_binary(model, cnn_inputs, SEQS, LABELS, MAX_LEN, epochs=60, seed=1)
    scores = predict_scores(model, cnn_inputs, SEQS, MAX_LEN)
    pos = scores[:len(AMP_FIXTURE)].mean()
    neg = scores[len(AMP_FIXTURE):].mean()
    assert pos > neg  # after training, positives outrank negatives on average
