import os.path as osp
import os
import random
import torch
import torch.nn.functional as F
import torch.nn as nn
import csv
from tqdm import tqdm
import numpy as np
import copy
import dgl
from torch_geometric.nn import GATConv
import argparse
from sklearn.metrics import (roc_auc_score, confusion_matrix, precision_recall_curve, 
                             average_precision_score, accuracy_score, f1_score, recall_score, precision_score)
from sklearn import preprocessing
import matplotlib.pyplot as plt
from utils.augmentation import shuffle_node_features
from torch_scatter import scatter_mean
from utils.gnn_explainer import GNNExplainer

import openpyxl
from networkx.algorithms import is_isomorphic

class AvgReadout(nn.Module):
    def __init__(self):
        super(AvgReadout, self).__init__()

    def forward(self, edge_index, seq, msk=None):
        row, col = edge_index

        c = scatter_mean(seq[row], col, dim=0, dim_size=seq.size(0))

        if msk is not None:
            msk = msk.unsqueeze(-1)
            c = c * msk
            c = c / msk.sum(dim=0).clamp(min=1e-6)
        return c

class Discriminator(nn.Module):
    def __init__(self, n_h, hidden_dim=32):
        super(Discriminator, self).__init__()

        self.bilinear1 = nn.Bilinear(n_h, n_h, hidden_dim)
        self.relu1 = nn.ReLU()

        self.fc1 = nn.Linear(hidden_dim, hidden_dim // 2)
        self.relu2 = nn.ReLU()

        self.fc2 = nn.Linear(hidden_dim // 2, 1)

        for m in self.modules():
            self.weights_init(m)

    def weights_init(self, m):
        if isinstance(m, nn.Bilinear) or isinstance(m, nn.Linear):
            torch.nn.init.xavier_uniform_(m.weight.data)
            if m.bias is not None:
                m.bias.data.fill_(0.0)

    def forward(self, c, h_pl, h_mi, s_bias1=None, s_bias2=None):

        sc_1 = self.bilinear1(h_pl, c)  # [N, hidden_dim]
        sc_1 = self.relu1(sc_1)

        sc_2 = self.bilinear1(h_mi, c)  # [N, hidden_dim]
        sc_2 = self.relu1(sc_2)

        sc_1 = self.fc1(sc_1)  # [N, hidden_dim // 2]
        sc_1 = self.relu2(sc_1)

        sc_2 = self.fc1(sc_2)  # [N, hidden_dim // 2]
        sc_2 = self.relu2(sc_2)

        sc_1 = self.fc2(sc_1).squeeze(1)  # [N]
        sc_2 = self.fc2(sc_2).squeeze(1)  # [N]

        if s_bias1 is not None:
            sc_1 += s_bias1
        if s_bias2 is not None:
            sc_2 += s_bias2

        logits = torch.cat((sc_1, sc_2), 0)  # [2N]

        return logits

class DGI(nn.Module):
    def __init__(self, n_in, n_h):

        super(DGI, self).__init__()
        self.fc_layer = nn.Linear(n_in, n_h)
        self.gat1 = GATConv(
            in_channels=n_in, 
            out_channels=n_h // 2, 
            heads=4, 
            activation=F.elu, 
            concat=True,
            add_self_loops=True, 
            negative_slope=0.2
        )
        self.gat2 = GATConv(
            in_channels=(n_h // 2) ,
            out_channels=n_h, 
            heads=4, 
            activation=F.elu, 
            concat=True, 
            add_self_loops=True, 
            negative_slope=0.2
        )

        self.read = AvgReadout()
        self.sigm = nn.Sigmoid()
        
        self.disc = Discriminator(40)

    def forward(self, x1, edge_index1, x2, edge_index2, samp_bias1=None, samp_bias2=None):

        h_fc_g1 = self.fc_layer(x1)  # [N, n_h]

        h_gat1_g1 = self.gat1(x1, edge_index1)
        if self.gat1.concat:
            h_gat1_g1 = h_gat1_g1.view(-1, self.gat1.heads, h_gat1_g1.size(-1) // self.gat1.heads).mean(dim=1) 

        h_gat2_g1 = self.gat2(h_gat1_g1, edge_index1)
        if self.gat2.concat:
            h_gat2_g1 = h_gat2_g1.view(-1, self.gat2.heads, h_gat2_g1.size(-1) // self.gat2.heads).mean(dim=1) 

        h1 = torch.cat([h_fc_g1, h_gat1_g1, h_gat2_g1], dim=1)

        c = self.read(edge_index1, h1)
        c = self.sigm(c)

        h_fc_g2 = self.fc_layer(x2) 

        h_gat1_g2 = self.gat1(x2, edge_index2)
        if self.gat1.concat:
            h_gat1_g2 = h_gat1_g2.view(-1, self.gat1.heads, h_gat1_g2.size(-1) // self.gat1.heads).mean(dim=1)


        h_gat2_g2 = self.gat2(h_gat1_g2, edge_index2) 
        if self.gat2.concat:
            h_gat2_g2 = h_gat2_g2.view(-1, self.gat2.heads, h_gat2_g2.size(-1) // self.gat2.heads).mean(dim=1) 



        h2 = torch.cat([h_fc_g2, h_gat1_g2, h_gat2_g2], dim=1) 

        logits = self.disc(c, h1, h2, samp_bias1, samp_bias2) 

        return logits

    def embed(self, x, edge_index, msk=None):

        h_fc = self.fc_layer(x)

        h_gat1 = self.gat1(x, edge_index) 
        if self.gat1.concat:
            h_gat1 = h_gat1.view(-1, self.gat1.heads, h_gat1.size(-1) // self.gat1.heads).mean(dim=1) 

        h_gat2 = self.gat2(h_gat1, edge_index) 
        if self.gat2.concat:
            h_gat2 = h_gat2.view(-1, self.gat2.heads, h_gat2.size(-1) // self.gat2.heads).mean(dim=1)

        h1 = torch.cat([h_fc, h_gat1, h_gat2], dim=1)

        c = self.read(edge_index, h1, msk=msk) 
        c = self.sigm(c)

        return h1.detach(), c.detach()

class LogReg(nn.Module):
    def __init__(self, ft_in, hidden_dim, nb_classes=1):
        super(LogReg, self).__init__()
        self.fc1 = nn.Linear(ft_in, hidden_dim)      
        self.relu = nn.ReLU()                        
        self.fc2 = nn.Linear(hidden_dim, nb_classes) 

        for m in self.modules():
            self.weights_init(m)

    def weights_init(self, m):
        if isinstance(m, nn.Linear):
            torch.nn.init.xavier_uniform_(m.weight.data)
            if m.bias is not None:
                m.bias.data.fill_(0.0)

    def forward(self, x):
        out = self.fc1(x)
        out = self.relu(out)
        out = self.fc2(out)
        return out

# @torch.no_grad()
def test(g_test_dgl, model, premod, device, epoch, test_chr):
    premod.eval()
    ys, preds = torch.tensor([]).to(device), torch.tensor([]).to(device)

    g_test_dgl = g_test_dgl.to(device)
    src1, dst1 = g_test_dgl.edges()
    edge_index = torch.stack([src1, dst1], dim=0).long().to(device)
    x = g_test_dgl.ndata['feat']
    y = g_test_dgl.ndata['label'].to(device)

    explainer = GNNExplainer(model, premod, epochs=30, lr=0.003, num_hops=2)

    all_edge_top5 = -np.ones((x.shape[0], 10))
    connected_edge_top5 = -np.ones((x.shape[0], 10))

    all_edge = [['1' for _ in range(14)] for _ in range(x.shape[0] + 10)]
    connected_edge = [['1' for _ in range(14)] for _ in range(x.shape[0] + 10)]

    all_edge[0][0] = 'node_index'
    all_edge[0][1] = 'edges'
    all_edge[0][11] = '1-order'
    all_edge[0][12] = '2-order'

    connected_edge[0][0] = 'node_index'
    for i in range(1, 11):
        connected_edge[0][i] = 'Edge'
    connected_edge[0][11] = '1-order'
    connected_edge[0][12] = '2-order'
    connected_edge[0][13] = 'Motif Mode'

    G_list = []
    topology_idx = []
    topology_num = []

    # print(x.shape[0])
    return 

    # for node_idx in range(0, x.shape[0]):
    for node_idx in range(0, x.shape[0]):
        if node_idx not in edge_index[0] and node_idx not in edge_index[1]:
            continue 

        node_feat_mask, edge_mask = explainer.explain_node(node_idx, x, edge_index)
        G, top5_G, connected_edge_list, connected_edge_importance = explainer.visualize_subgraph(node_idx, edge_index, edge_mask, y=y)
        edge_repeat = [node_idx]
        connected_edge_repeat = [node_idx]

        for i in range(len(connected_edge_list)):
            for j in range(len(connected_edge_list)):
                if (connected_edge_list[i][0], connected_edge_list[i][1]) == (connected_edge_list[j][1], connected_edge_list[j][0]):
                    if connected_edge_importance[i] > connected_edge_importance[j]:
                        connected_edge_importance[j] = 0
                    else:
                        connected_edge_importance[i] = 0

        if len(connected_edge_importance) != 0:
            connected_edge_importance[0] = connected_edge_importance[0][0]

        connected_edge_importance = torch.tensor(connected_edge_importance)

        if len(list(G.nodes)) == 1:
            all_edge_top5[node_idx][0] = node_idx
            continue

        sorted_edge_mask = torch.sort(edge_mask, descending=True)
        sorted_connected_edge_mask = torch.sort(connected_edge_importance, descending=True)
        all_edge[node_idx + 1][0] = str(node_idx + 1)
        connected_edge[node_idx + 1][0] = str(node_idx + 1)
        hop_1 = ''
        hop_2_list = []
        hop_2 = ''

        all_edge[node_idx + 1][11] = hop_1
        all_edge[node_idx + 1][12] = hop_2

        hop_1 = ''
        hop_1_list = []
        hop_2_list = []
        hop_2 = ''

        for i in range(0, min(10, len(connected_edge_list))):

            index = sorted_connected_edge_mask.indices[i]
            edge_1, edge_2 = connected_edge_list[index][0], connected_edge_list[index][1]

            importance = sorted_connected_edge_mask.values[i].item()
            if importance == 0:
                break

            all_str_edge = "(" + str(min(edge_1, edge_2).item() + 1) + ","
            all_str_edge = all_str_edge + str(max(edge_1, edge_2).item() + 1) + ") " + str(format(importance, '.5f'))

            connected_edge[node_idx + 1][i + 1] = all_str_edge
            if edge_1.item() == node_idx and edge_2.item() not in hop_1_list:
                hop_1_list.append(edge_2.item())
                hop_1 = hop_1 + ' ' + str(edge_2.item() + 1)
            elif edge_2.item() == node_idx and edge_1.item() not in hop_1_list:
                hop_1_list.append(edge_1.item())
                hop_1 = hop_1 + ' ' + str(edge_1.item() + 1)
            else:
                if edge_1.item() not in hop_2_list:
                    hop_2_list.append(edge_1.item())
                    hop_2 = hop_2 + ' ' + str(edge_1.item() + 1)
                elif edge_2.item() not in hop_2_list:
                    hop_2_list.append(edge_2.item())
                    hop_2 = hop_2 + ' ' + str(edge_2.item() + 1)


        connected_edge[node_idx + 1][11] = hop_1
        connected_edge[node_idx + 1][12] = hop_2

        new_topology = True

        if G_list == []:
            G_list.append(top5_G)
            topology_idx.append(node_idx)
            motif_num = 1
            topology_num.append(1)
            plt.title("motif = " + str(len(topology_num)))
            # plt.show()
            plt.savefig('explainer/NHEK_motif_' + str(test_chr) + 'motif_'  + str(len(topology_num))+'.png')
            plt.close()
            continue

        for i in range(len(G_list)):
            if is_isomorphic(G_list[i], top5_G) == True:
                topology_num[i] += 1
                motif_num = i + 1
                new_topology = False
                plt.close()
                break

        if new_topology == True:
            G_list.append(top5_G)
            topology_idx.append(node_idx)
            topology_num.append(1)
            motif_num = len(topology_num)
            plt.title("motif = " + str(len(topology_num)))
            # plt.show()
            plt.savefig('explainer/NHEK_motif_' + str(test_chr) + '/motif_'  + str(len(topology_num))+'.png')
            plt.close()

        connected_edge[node_idx + 1][13] = str(motif_num)

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = 'connected_edge'
    for r in range(len(connected_edge)):
        for c in range(len(connected_edge[0])):
            ws.cell(r + 1, c + 1).value = connected_edge[r][c]
    wb.save('explainer/NHEK_chr_' + str(test_chr) +'_connected_edge.xlsx')

    print('zero_topology_num:%s' % (topology_num))
    print('zero_topology_idx:%s' % (topology_idx))


parser = argparse.ArgumentParser(description='GAT Binary Classification Training with Dual Graphs')
parser.add_argument('--device', type=int, default=0, help='GPU device')
parser.add_argument('--learning_rate', type=float, default=0.01, help='Learning rate')
parser.add_argument('--epoch', type=int, default=300, help='Number of training epochs')

args = parser.parse_args()

device = torch.device('cuda:' + str(args.device) if torch.cuda.is_available() else 'cpu')
# device = torch.device('cpu')

for test_chr in range(1, 23):
    os.makedirs(f'explainer/NHEK_motif_{test_chr}', exist_ok=True)
    print(f'----------------------------------------')
    print(f'   Test_chr:{test_chr}')
    print(f'----------------------------------------')

    data_dgl = []
    # data_sdgl = []

    for i in range(1, 23):
        chr_str = f"chr{i}"

        g_dgl = dgl.load_graphs(f"preprocess/DGL/{chr_str}.dgl")[0][0]
        # g_sdgl = dgl.load_graphs(f"preprocess/SDGL/{chr_str}.dgl")[0][0]
        
        g_feats = np.loadtxt('preprocess/feats/chr' + str(i) + '_features.txt', delimiter='\t')
        g_dgl.ndata['feat'] = torch.from_numpy(g_feats).float()
        # g_sdgl.ndata['feat'] = torch.from_numpy(g_feats).float()

        data_dgl.append(g_dgl)
        # data_sdgl.append(g_sdgl)

    if test_chr == 22:
        valid_chr = 1
        g_test_dgl = data_dgl[21]
        # g_test_sdgl = data_sdgl[21]
        data_dgl.pop(21)
        # data_sdgl.pop(21)
        g_valid_dgl = data_dgl[0]
        # g_valid_sdgl = data_sdgl[0]
        data_dgl.pop(0)
        # data_sdgl.pop(0)
    else:
        valid_chr = test_chr + 1
        g_test_dgl = data_dgl[test_chr - 1]
        # g_test_sdgl = data_sdgl[test_chr - 1]
        g_valid_dgl = data_dgl[valid_chr - 1]
        # g_valid_sdgl = data_sdgl[valid_chr - 1]
        data_dgl.pop(test_chr - 1)
        # data_sdgl.pop(test_chr - 1)
        data_dgl.pop(test_chr - 1)
        # data_sdgl.pop(test_chr - 1)

    train_node_feats = torch.zeros(0, 5)
    train_edge_feats = torch.zeros(0, 1)

    for graph_dgl in data_dgl:
        train_node_feats = torch.cat([train_node_feats, graph_dgl.ndata['feat']], 0)
        train_edge_feats = torch.cat([train_edge_feats, graph_dgl.edata['edge_feature'].view(graph_dgl.edata['edge_feature'].shape[0], 1)], 0)

    density_scaler = preprocessing.StandardScaler().fit(train_node_feats)
    edge_scaler = preprocessing.StandardScaler().fit(train_edge_feats)

    for graph_dgl in data_dgl:
        graph_dgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(graph_dgl.ndata['feat'].cpu().numpy())).float()
        # graph_sdgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(graph_sdgl.ndata['feat'].cpu().numpy())).float()
        
        graph_dgl.edata['edge_feature'] = graph_dgl.edata['edge_feature'].view(graph_dgl.edata['edge_feature'].shape[0], 1)
        graph_dgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(graph_dgl.edata['edge_feature'])).float()
        # graph_sdgl.edata['edge_feature'] = graph_sdgl.edata['edge_feature'].view(graph_sdgl.edata['edge_feature'].shape[0], 1)
        # graph_sdgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(graph_sdgl.edata['edge_feature'])).float()

    g_valid_dgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(g_valid_dgl.ndata['feat'].cpu().numpy())).float()
    # g_valid_sdgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(g_valid_sdgl.ndata['feat'].cpu().numpy())).float()

    g_valid_dgl.edata['edge_feature'] = g_valid_dgl.edata['edge_feature'].view(g_valid_dgl.edata['edge_feature'].shape[0], 1)
    g_valid_dgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(g_valid_dgl.edata['edge_feature'])).float()
    # g_valid_sdgl.edata['edge_feature'] = g_valid_sdgl.edata['edge_feature'].view(g_valid_sdgl.edata['edge_feature'].shape[0], 1)
    # g_valid_sdgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(g_valid_sdgl.edata['edge_feature'])).float()

    g_test_dgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(g_test_dgl.ndata['feat'].cpu().numpy())).float()
    # g_test_sdgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(g_test_sdgl.ndata['feat'].cpu().numpy())).float()

    g_test_dgl.edata['edge_feature'] = g_test_dgl.edata['edge_feature'].view(g_test_dgl.edata['edge_feature'].shape[0], 1)
    g_test_dgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(g_test_dgl.edata['edge_feature'])).float()
    # g_test_sdgl.edata['edge_feature'] = g_test_sdgl.edata['edge_feature'].view(g_test_sdgl.edata['edge_feature'].shape[0], 1)
    # g_test_sdgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(g_test_sdgl.edata['edge_feature'])).float()

    model = DGI(5, 16).to(device)
    model.load_state_dict(torch.load('result/model/best_dgi2.pth'))

    premod = LogReg(40, 32).to(device)
    premod.load_state_dict(torch.load('result/model/best_premod2.pth'))

    test(g_test_dgl, model, premod, device, args.epoch, test_chr)



