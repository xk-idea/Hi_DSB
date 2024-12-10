import dgl
import torch
import numpy as np
import copy

def shuffle_node_features(g, feature_key='feat'):

    g_shuffled = copy.deepcopy(g)
    

    if feature_key not in g_shuffled.ndata:
        raise KeyError(f"Key '{feature_key}' not exist in the node data of the graph")
    
    features = g_shuffled.ndata[feature_key]

    if features.dim() == 2:
        num_nodes = features.size(0)
        perm = torch.randperm(num_nodes)
        g_shuffled.ndata[feature_key] = features[perm]
    else:
        raise ValueError("The supported feature dimension is 2, but the received feature dimension is {features.dim()}。")
    
    return g_shuffled

