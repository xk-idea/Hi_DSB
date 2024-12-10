import os.path as osp
import os
import random
import torch
import torch.nn.functional as F
import torch.nn as nn
import csv
import time
from tqdm import tqdm
import numpy as np
import copy
from sklearn import svm
import dgl
from dgl.nn import GATConv  
import argparse
from sklearn.metrics import (roc_auc_score, confusion_matrix, precision_recall_curve, 
                             average_precision_score, accuracy_score, f1_score, recall_score, precision_score)
from sklearn import preprocessing
from utils.augmentation import shuffle_node_features

class AvgReadout(nn.Module):
    def __init__(self):
        super(AvgReadout, self).__init__()

    def forward(self, graph, seq, msk=None):

        with graph.local_scope():
            graph.ndata['h'] = seq

            graph.update_all(message_func=dgl.function.copy_u('h', 'm'),
                            reduce_func=dgl.function.mean('m', 'c'))
            c = graph.ndata['c']  

            if msk is not None:
                msk = msk.unsqueeze(-1)  
                c = c * msk
                c = c / msk.sum(dim=0)  
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


        sc_1 = self.bilinear1(h_pl, c)
        sc_1 = self.relu1(sc_1)        

        sc_2 = self.bilinear1(h_mi, c)
        sc_2 = self.relu1(sc_2)

        sc_1 = self.fc1(sc_1)
        sc_1 = self.relu2(sc_1)
        sc_2 = self.fc1(sc_2)
        sc_2 = self.relu2(sc_2)

        sc_1 = self.fc2(sc_1).squeeze(1)
        sc_2 = self.fc2(sc_2).squeeze(1)
        
        if s_bias1 is not None:
            sc_1 += s_bias1
        if s_bias2 is not None:
            sc_2 += s_bias2

        logits = torch.cat((sc_1, sc_2), 0)
        return logits

class DGI(nn.Module):
    def __init__(self, n_in, n_h):
        super(DGI, self).__init__()

        self.fc_layer = nn.Linear(n_in, n_h)
        
        self.gat1 = GATConv(n_in, n_h // 2, num_heads=4, activation=F.elu, allow_zero_in_degree=True)
        self.gat2 = GATConv(n_h // 2, n_h, num_heads=4, activation=F.elu, allow_zero_in_degree=True)
        
        self.read = AvgReadout()
        self.sigm = nn.Sigmoid()
        
        # self.disc = Discriminator(40)
        self.disc = Discriminator(n_h)
        
    def forward(self, g1, g2, samp_bias1, samp_bias2):

        h_fc_g1 = self.fc_layer(g1.ndata['feat'])  # [N, n_h]
        h_gat1_g1 = self.gat1(g1, g1.ndata['feat']).mean(dim=1)  # [N, n_h]
        h_gat2_g1 = self.gat2(g1, h_gat1_g1).mean(dim=1)  # [N, n_h]
        h1 = torch.cat([h_fc_g1, h_gat1_g1, h_gat2_g1], dim=1)  # [N, 3 *n_h]

        c = self.read(g1, h1)  # [N, 3 *n_h]
        c = self.sigm(c)

        h_fc_g2 = self.fc_layer(g2.ndata['feat'])  # [N, n_h]
        h_gat1_g2 = self.gat1(g2, g2.ndata['feat']).mean(dim=1)  # [N, n_h]
        h_gat2_g2 = self.gat2(g2, h_gat1_g2).mean(dim=1)  # [N, n_h]
        h2 = torch.cat([h_fc_g2, h_gat1_g2, h_gat2_g2], dim=1)  # [N, 3 *n_h]

        logits = self.disc(c, h1, h2, samp_bias1, samp_bias2)

        return logits

    def embed(self, g1, msk):
        h_fc = self.fc_layer(g1.ndata['feat'])
        h_gat1 = self.gat1(g1, g1.ndata['feat']).mean(dim=1)
        h_gat2 = self.gat2(g1, h_gat1).mean(dim=1)
        h1 = torch.cat([h_fc, h_gat1, h_gat2], dim=1)
        c = self.read(g1, h1)
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

def train(epoch, train_list, model, premod, optimizer, device):
    premod.train()

    length = len(train_list)
    pbar = tqdm(total=length)
    pbar.set_description(f'Epoch {epoch:02d}')

    total_loss = 0
    ys, preds = torch.tensor([]).to(device), torch.tensor([]).to(device)

    loss_op = torch.nn.BCEWithLogitsLoss()

    for i, (g_train_dgl, g_train_sdgl) in enumerate(train_list):
        g_train_dgl = g_train_dgl.to(device)
        g_train_sdgl = g_train_sdgl.to(device)


        x1, _ = model.embed(g_train_dgl, None)
        x2, _ = model.embed(g_train_sdgl, None)
        x = torch.cat([x1, x2], dim=1).to(device)
        y = g_train_dgl.ndata['label'].to(device)

        optimizer.zero_grad()

        logits = premod(x)
        loss = loss_op(logits.squeeze(), y.squeeze().float())

        total_loss += loss.item()
        loss.backward()
        optimizer.step()

        ys = torch.cat([ys, y])
        preds = torch.cat([preds, logits])

        pbar.update(1)

    pbar.close()
    ys, preds = ys.squeeze(), preds.squeeze()
    ys_np, preds_np = ys.cpu().detach().numpy(), preds.cpu().detach().numpy()


    binary_preds = (preds_np >= 0).astype(int)
    acc = accuracy_score(ys_np, binary_preds)
    f1 = f1_score(ys_np, binary_preds)
    recall = recall_score(ys_np, binary_preds)
    precision = precision_score(ys_np, binary_preds)

    return (total_loss / length, 
            roc_auc_score(ys_np, preds_np), 
            average_precision_score(ys_np, preds_np),
            acc, f1, recall, precision, 
            ys_np, preds_np)

@torch.no_grad()
def validation(g_valid_dgl, g_valid_sdgl, model, premod, device):
    premod.eval()
    
    loss_op = torch.nn.BCEWithLogitsLoss()

    ys, preds = torch.tensor([]).to(device), torch.tensor([]).to(device)

    g_valid_dgl = g_valid_dgl.to(device)
    g_valid_sdgl = g_valid_sdgl.to(device)

    x1, _ = model.embed(g_valid_dgl, None)
    x2, _ = model.embed(g_valid_sdgl, None)
    x = torch.cat([x1, x2], dim=1).to(device)
    y = g_valid_dgl.ndata['label'].to(device)

    logits = premod(x)
    loss = loss_op(logits.squeeze(), y.squeeze().float())

    ys = torch.cat([ys, y])
    preds = torch.cat([preds, logits])

    ys, preds = ys.squeeze(), preds.squeeze()
    ys_np, preds_np = ys.cpu().detach().numpy(), preds.cpu().detach().numpy()


    binary_preds = (preds_np >= 0).astype(int)
    acc = accuracy_score(ys_np, binary_preds)
    f1 = f1_score(ys_np, binary_preds)
    recall = recall_score(ys_np, binary_preds)
    precision = precision_score(ys_np, binary_preds)

    return (roc_auc_score(ys_np, preds_np), 
            average_precision_score(ys_np, preds_np), 
            acc, f1, recall, precision, 
            ys_np, preds_np, loss.item())

@torch.no_grad()
def test(g_test_dgl, g_test_sdgl, model, premod, device):
    premod.eval()
    ys, preds = torch.tensor([]).to(device), torch.tensor([]).to(device)

    g_test_dgl = g_test_dgl.to(device)
    g_test_sdgl = g_test_sdgl.to(device)

    x1, _ = model.embed(g_test_dgl, None)
    x2, _ = model.embed(g_test_sdgl, None)
    x = torch.cat([x1, x2], dim=1).to(device)

    y = g_test_dgl.ndata['label'].to(device)

    logits = premod(x)

    ys = torch.cat([ys, y])
    preds = torch.cat([preds, logits])

    ys, preds = ys.squeeze(), preds.squeeze()
    ys_np, preds_np = ys.cpu().detach().numpy(), preds.cpu().detach().numpy()


    binary_preds = (preds_np >= 0).astype(int)
    acc = accuracy_score(ys_np, binary_preds)
    f1 = f1_score(ys_np, binary_preds)
    recall = recall_score(ys_np, binary_preds)
    precision = precision_score(ys_np, binary_preds)

    return (roc_auc_score(ys_np, preds_np), 
            average_precision_score(ys_np, preds_np), 
            acc, f1, recall, precision, 
            ys_np, preds_np)

parser = argparse.ArgumentParser(description='GAT Binary Classification Training with Dual Graphs')
parser.add_argument('--device', type=int, default=7, help='GPU device')
parser.add_argument('--learning_rate', type=float, default=0.001, help='Learning rate')
parser.add_argument('--epoch', type=int, default=300, help='Number of training epochs')

args = parser.parse_args()

device = torch.device('cuda:' + str(args.device) if torch.cuda.is_available() else 'cpu')
# device = torch.device('cpu')


for test_chr in range(1, 23):
    print(f'----------------------------------------')
    print(f'   Test_chr:{test_chr}')
    print(f'----------------------------------------')

    train_list = []
    data_dgl = []
    data_sdgl = []

    for i in range(1, 23):
        chr_str = f"chr{i}"


        g_dgl = dgl.load_graphs(f"preprocess/DGL/{chr_str}.dgl")[0][0]

        g_sdgl = dgl.load_graphs(f"preprocess/SDGL/{chr_str}.dgl")[0][0]
        

        g_feats = np.loadtxt('preprocess/4feats/chr' + str(i) + '_features.txt', delimiter='\t')


        g_dgl.ndata['feat'] = torch.from_numpy(g_feats).float()
        g_sdgl.ndata['feat'] = torch.from_numpy(g_feats).float()

        data_dgl.append(g_dgl)
        data_sdgl.append(g_sdgl)


    if test_chr == 22:
        valid_chr = 1
        g_test_dgl = data_dgl[21]
        g_test_sdgl = data_sdgl[21]
        data_dgl.pop(21)
        data_sdgl.pop(21)
        g_valid_dgl = data_dgl[0]
        g_valid_sdgl = data_sdgl[0]
        data_dgl.pop(0)
        data_sdgl.pop(0)
    else:
        valid_chr = test_chr + 1
        g_test_dgl = data_dgl[test_chr - 1]
        g_test_sdgl = data_sdgl[test_chr - 1]
        g_valid_dgl = data_dgl[valid_chr - 1]
        g_valid_sdgl = data_sdgl[valid_chr - 1]
        data_dgl.pop(test_chr - 1)
        data_sdgl.pop(test_chr - 1)
        data_dgl.pop(test_chr - 1)
        data_sdgl.pop(test_chr - 1)


    train_node_feats = torch.zeros(0, 4)
    train_edge_feats = torch.zeros(0, 1)
    
    for graph_dgl, graph_sdgl in zip(data_dgl, data_sdgl):
        train_node_feats = torch.cat([train_node_feats, graph_dgl.ndata['feat']], 0)
        train_edge_feats = torch.cat([train_edge_feats, graph_dgl.edata['edge_feature'].view(graph_dgl.edata['edge_feature'].shape[0], 1)], 0)

    density_scaler = preprocessing.StandardScaler().fit(train_node_feats)
    edge_scaler = preprocessing.StandardScaler().fit(train_edge_feats)

    for graph_dgl, graph_sdgl in zip(data_dgl, data_sdgl):
        graph_dgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(graph_dgl.ndata['feat'].cpu().numpy())).float()
        graph_sdgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(graph_sdgl.ndata['feat'].cpu().numpy())).float()
        
        graph_dgl.edata['edge_feature'] = graph_dgl.edata['edge_feature'].view(graph_dgl.edata['edge_feature'].shape[0], 1)
        graph_dgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(graph_dgl.edata['edge_feature'])).float()
        graph_sdgl.edata['edge_feature'] = graph_sdgl.edata['edge_feature'].view(graph_sdgl.edata['edge_feature'].shape[0], 1)
        graph_sdgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(graph_sdgl.edata['edge_feature'])).float()

    g_valid_dgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(g_valid_dgl.ndata['feat'].cpu().numpy())).float()
    g_valid_sdgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(g_valid_sdgl.ndata['feat'].cpu().numpy())).float()
    
    g_valid_dgl.edata['edge_feature'] = g_valid_dgl.edata['edge_feature'].view(g_valid_dgl.edata['edge_feature'].shape[0], 1)
    g_valid_dgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(g_valid_dgl.edata['edge_feature'])).float()
    g_valid_sdgl.edata['edge_feature'] = g_valid_sdgl.edata['edge_feature'].view(g_valid_sdgl.edata['edge_feature'].shape[0], 1)
    g_valid_sdgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(g_valid_sdgl.edata['edge_feature'])).float()

    g_test_dgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(g_test_dgl.ndata['feat'].cpu().numpy())).float()
    g_test_sdgl.ndata['feat'] = torch.from_numpy(density_scaler.transform(g_test_sdgl.ndata['feat'].cpu().numpy())).float()
    
    g_test_dgl.edata['edge_feature'] = g_test_dgl.edata['edge_feature'].view(g_test_dgl.edata['edge_feature'].shape[0], 1)
    g_test_dgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(g_test_dgl.edata['edge_feature'])).float()
    g_test_sdgl.edata['edge_feature'] = g_test_sdgl.edata['edge_feature'].view(g_test_sdgl.edata['edge_feature'].shape[0], 1)
    g_test_sdgl.edata['feat'] = torch.from_numpy(edge_scaler.transform(g_test_sdgl.edata['edge_feature'])).float()


    train_list = list(zip(data_dgl, data_sdgl))

    model = DGI(4, 16).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.learning_rate)
    AL = nn.BCEWithLogitsLoss()

    best_loss = inf

    for epoch1 in range(20):
        model.train()
        for i, (g_train_dgl, g_train_sdgl) in enumerate(train_list):
            g1 = g_train_dgl.to(device)
            g2 = g_train_sdgl.to(device)
            g1_aug=shuffle_node_features(g1)
            g2_aug=shuffle_node_features(g2)

            lbl_1 = torch.ones(g1.ndata['feat'].size(0), 1)
            lbl_2 = torch.zeros(g1.ndata['feat'].size(0), 1)
            lbl = torch.cat((lbl_1, lbl_2), 0).squeeze().to(device)
            optimizer.zero_grad()
            logits = model(g1, g1_aug, None, None)
            loss = AL(logits, lbl)
            if loss < best_loss:
                best_loss = loss
                torch.save(model.state_dict(), 'result/model/best_dgi.pkl')
            loss.backward()
            optimizer.step()
            print(f'{epoch1+1}.Loss:{loss}')


    model.load_state_dict(torch.load('result/model/best_dgi.pkl'))
    premod = LogReg(80, 32).to(device)
    optimizer_premod = torch.optim.Adam(premod.parameters(), lr=args.learning_rate, weight_decay=5e-4)


    best_loss_val = float('inf')
    best_epoch = -1
    best_valid_auc = -1
    best_valid_prauc = -1
    best_valid_acc = -1
    best_valid_f1 = -1
    best_valid_recall = -1
    best_valid_precision = -1
    best_test_auc = -1
    best_test_prauc = -1
    best_test_acc = -1
    best_test_f1 = -1
    best_test_recall = -1
    best_test_precision = -1
    es = 0  

    for epoch in range(1, args.epoch + 1):
        # shuffle
        random.shuffle(train_list)

        result_train = train(
            epoch, train_list, model, premod, optimizer_premod, device
        )
        (loss, train_auc, train_prauc, train_acc, train_f1, 
            train_recall, train_precision, 
            train_y_train, train_y_pred) = result_train

        result_valid = validation(
            g_valid_dgl, g_valid_sdgl, model, premod, device
        )
        (valid_auc, valid_prauc, valid_acc, valid_f1, 
            valid_recall, valid_precision, 
            valid_y_test, valid_y_pred, valid_loss) = result_valid

        result_test = test(
            g_test_dgl, g_test_sdgl, model, premod, device
        )
        (test_auc, test_prauc, test_acc, test_f1, 
            test_recall, test_precision, 
            test_y_test, test_y_pred) = result_test

        print('\n')
        result_str = (f'test_chr: {test_chr:02d}: '
                        f'Epoch: {epoch:02d}, Train_Loss: {loss:.4f}, '
                        f'Train_AUC: {train_auc:.4f}, Train_PRC: {train_prauc:.4f}, '
                        f'Train_ACC: {train_acc:.4f}, Train_F1: {train_f1:.4f}, '
                        f'Train_Recall: {train_recall:.4f}, Train_Precision: {train_precision:.4f};   '
                        f'Valid_Loss: {valid_loss:.4f}, Valid_AUC: {valid_auc:.4f}, '
                        f'Valid_PRC: {valid_prauc:.4f}, Valid_ACC: {valid_acc:.4f}, '
                        f'Valid_F1: {valid_f1:.4f}, Valid_Recall: {valid_recall:.4f}, '
                        f'Valid_Precision: {valid_precision:.4f};   '
                        f'Test_AUC: {test_auc:.4f}, Test_PRC: {test_prauc:.4f}, '
                        f'Test_ACC: {test_acc:.4f}, Test_F1: {test_f1:.4f}, '
                        f'Test_Recall: {test_recall:.4f}, Test_Precision: {test_precision:.4f}')
        print(result_str)

        if valid_loss < best_loss_val:
            best_loss_val = valid_loss

            best_epoch = epoch
            best_valid_auc = valid_auc
            best_valid_prauc = valid_prauc
            best_valid_acc = valid_acc
            best_valid_f1 = valid_f1
            best_valid_recall = valid_recall
            best_valid_precision = valid_precision
            best_test_auc = test_auc
            best_test_prauc = test_prauc
            best_test_acc = test_acc
            best_test_f1 = test_f1
            best_test_recall = test_recall
            best_test_precision = test_precision
            checkpoint = copy.deepcopy(model)
            es = 0
        else:
            es += 1
            if es > 9:
                print("Early stopping triggered.")
                break

    best_str = (f'test_chr: {test_chr:02d}: '
                f'Best_Epoch:{best_epoch:02d}, Best_Valid_Loss:{best_loss_val:.4f}, '
                f'Best_Valid_AUC:{best_valid_auc:.4f}, Best_Valid_PRC:{best_valid_prauc:.4f}, '
                f'Best_Valid_ACC:{best_valid_acc:.4f}, Best_Valid_F1:{best_valid_f1:.4f}, '
                f'Best_Valid_Recall:{best_valid_recall:.4f}, Best_Valid_Precision:{best_valid_precision:.4f}; '
                f'Best_Test_AUC:{best_test_auc:.4f}, Best_Test_PRC:{best_test_prauc:.4f}, '
                f'Best_Test_ACC:{best_test_acc:.4f}, Best_Test_F1:{best_test_f1:.4f}, '
                f'Best_Test_Recall:{best_test_recall:.4f}, Best_Test_Precision:{best_test_precision:.4f}')
    print(best_str)

    train_type_name = '4feats_DGL'
    os.makedirs('result/precision/', exist_ok=True)
    with open(f'result/precision/{train_type_name}.txt', 'a') as file:
        file.write(best_str + '\n')