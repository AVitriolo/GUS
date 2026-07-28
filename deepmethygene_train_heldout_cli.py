from scipy.stats import pearsonr
import pandas as pd
import math
from tqdm import tqdm
import numpy as np
from sklearn.model_selection import train_test_split, KFold
import os, logging, pickle, random, torch
from sklearn.metrics import r2_score
import torch.nn as nn
import warnings
import gc
import argparse

warnings.filterwarnings("ignore")
SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed_all(42)

parser = argparse.ArgumentParser()
parser.add_argument("--input_path_dataset_train", required=True, type=str)
parser.add_argument("--input_path_dataset_test",  required=True, type=str)
parser.add_argument("--tx_id",       required=True, type=str)
parser.add_argument("--output_file", required=True, type=str)
args = parser.parse_args()

device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# --- Load pre-split matrices (same TSVs model_selection_split_kbest.py consumes) ---
# Split happens upstream (split_samples.py); do NOT re-split here.
TRAIN = pd.read_csv(args.input_path_dataset_train, header=0, sep="\t", index_col=0)
TEST  = pd.read_csv(args.input_path_dataset_test,  header=0, sep="\t", index_col=0)

TxID = args.tx_id
assert TxID in TRAIN.columns and TxID in TEST.columns, f"{TxID} missing from a matrix"
assert list(TRAIN.columns) == list(TEST.columns), "train/test columns differ"

feature_cols = [c for c in TRAIN.columns if c != TxID]
train_pool_df = TRAIN[feature_cols + [TxID]].sample(frac=1, random_state=1)
held_out_df   = TEST[feature_cols + [TxID]]

assert not (set(train_pool_df.index) & set(held_out_df.index)), "train/test sample overlap"

train_pool_x = train_pool_df[feature_cols].values
train_pool_y = train_pool_df[TxID].values
held_out_x   = held_out_df[feature_cols].values
held_out_y   = held_out_df[TxID].values

print(f"[{TxID}] train-pool {train_pool_x.shape[0]}x{train_pool_x.shape[1]} | held-out {held_out_x.shape[0]}x{held_out_x.shape[1]}")

class ResidualBlock(nn.Module):
    def __init__(self, channels):
        super(ResidualBlock, self).__init__()
        self.conv1 = nn.Conv1d(channels, channels, kernel_size=3, padding=1)
        self.conv2 = nn.Conv1d(channels, channels, kernel_size=3, padding=1)
        self.leaky_relu = nn.LeakyReLU(0.01)

    def forward(self, x):
        residual = x
        out = self.leaky_relu(self.conv1(x))
        out = self.conv2(out)
        out += residual
        out = self.leaky_relu(out)
        return out

class AdaptiveRegressionCNN(nn.Module):
    def __init__(self, input_size):
        super(AdaptiveRegressionCNN, self).__init__()
        out_channels_conv1 = max(1, min(64, input_size // 10))
        out_channels_conv2 = max(1, min(32, input_size // 20))
        self.conv1 = nn.Conv1d(in_channels=1, out_channels=out_channels_conv1, kernel_size=3, padding=1)
        self.resblock1 = ResidualBlock(out_channels_conv1)
        self.conv2 = nn.Conv1d(in_channels=out_channels_conv1, out_channels=out_channels_conv2, kernel_size=3, padding=1)
        self.resblock2 = ResidualBlock(out_channels_conv2)
        self.leaky_relu = nn.LeakyReLU(0.01)
        self._to_linear = None
        self._calculate_to_linear(input_size)
        self.fc1 = nn.Linear(self._to_linear, 512)
        self.fc2 = nn.Linear(512, 1)

    def _calculate_to_linear(self, L):
        x = torch.randn(1, 1, L)
        x = self.leaky_relu(self.conv1(x))
        x = self.resblock1(x)
        x = self.leaky_relu(self.conv2(x))
        x = self.resblock2(x)
        self._to_linear = x.numel() // x.size(0)

    def forward(self, x):
        x = self.leaky_relu(self.conv1(x))
        x = self.resblock1(x)
        x = self.leaky_relu(self.conv2(x))
        x = self.resblock2(x)
        x = x.view(x.size(0), -1)
        x = self.leaky_relu(self.fc1(x))
        x = self.fc2(x)
        return x


def train_one_model(model, x_train, y_train, x_val, y_val, patience, epochs, device):
    """Train with early stopping on (x_val, y_val); returns the best-epoch model."""
    model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.MSELoss()
    best_val_loss = float('inf')
    count = 0
    best_model_state = model.state_dict()

    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        y_pred = model(x_train)
        loss = criterion(y_pred, y_train)
        loss.backward()
        optimizer.step()

        model.eval()
        with torch.no_grad():
            y_pred_val = model(x_val)
            val_loss = criterion(y_pred_val, y_val)
            if val_loss.item() < best_val_loss:
                best_val_loss = val_loss.item()
                count = 0
                best_model_state = model.state_dict()
            else:
                count += 1
                if count >= patience:
                    break

    model.load_state_dict(best_model_state)
    return model


def cross_validation(Dmodel_ctor, x, y, patience, epochs, device, n_splits=5):
    """5-fold CV WITHIN whatever data it's given (here: train pool only)."""
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    r2_scores = []
    cv_r2_scores = []
    for train_idx, val_idx in kf.split(x):
        x_train, x_val = x[train_idx], x[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        x_train_t = torch.tensor(x_train, dtype=torch.float32).unsqueeze(1).to(device)
        y_train_t = torch.tensor(y_train, dtype=torch.float32).unsqueeze(1).to(device)
        x_val_t = torch.tensor(x_val, dtype=torch.float32).unsqueeze(1).to(device)
        y_val_t = torch.tensor(y_val, dtype=torch.float32).unsqueeze(1).to(device)

        model = Dmodel_ctor()
        model = train_one_model(model, x_train_t, y_train_t, x_val_t, y_val_t, patience, epochs, device)

        model.eval()
        with torch.no_grad():
            y_pred_val = model(x_val_t)
            last_pre = [i[0] for i in y_pred_val.cpu().numpy().tolist()]
            last_val = [i[0] for i in y_val_t.cpu().numpy().tolist()]
            r = pearsonr(last_pre, last_val)
            r2 = math.pow(r[0], 2)
            r2_scores.append(r2)
            cv_r2 = r2_score(last_val, last_pre)
            cv_r2_scores.append(cv_r2)
        del model, x_train_t, x_val_t, y_train_t, y_val_t, y_pred_val
        torch.cuda.empty_cache()
        gc.collect()
    return np.mean(r2_scores), np.mean(cv_r2_scores)


input_size = train_pool_x.shape[1]
torch.manual_seed(42)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(42)

def model_ctor():
    torch.manual_seed(42)
    return AdaptiveRegressionCNN(input_size=input_size)

# --- 1. 5-fold CV within the train pool only (same style as other DeepMethyGene runs, for comparability) ---
mean_r2, mean_cv_r2 = cross_validation(model_ctor, train_pool_x, train_pool_y, patience=90, epochs=700, device=device)

# --- 2. Final model: train on train pool with an internal (train-pool-only) validation split for early stopping ---
inner_train_x, inner_val_x, inner_train_y, inner_val_y = train_test_split(
    train_pool_x, train_pool_y, test_size=0.15, random_state=1
)
inner_train_x_t = torch.tensor(inner_train_x, dtype=torch.float32).unsqueeze(1).to(device)
inner_train_y_t = torch.tensor(inner_train_y, dtype=torch.float32).unsqueeze(1).to(device)
inner_val_x_t = torch.tensor(inner_val_x, dtype=torch.float32).unsqueeze(1).to(device)
inner_val_y_t = torch.tensor(inner_val_y, dtype=torch.float32).unsqueeze(1).to(device)

final_model = model_ctor()
final_model = train_one_model(final_model, inner_train_x_t, inner_train_y_t, inner_val_x_t, inner_val_y_t,
                                patience=90, epochs=700, device=device)

# --- 3. Evaluate the final model on the FULL train pool (train_r2, mirrors GUS's train_r2) ---
train_pool_x_t = torch.tensor(train_pool_x, dtype=torch.float32).unsqueeze(1).to(device)
train_pool_y_t = torch.tensor(train_pool_y, dtype=torch.float32).unsqueeze(1).to(device)
final_model.eval()
with torch.no_grad():
    y_pred_train = final_model(train_pool_x_t)
    train_pred_list = [v[0] for v in y_pred_train.cpu().numpy().tolist()]
    train_true_list = [v[0] for v in train_pool_y_t.cpu().numpy().tolist()]
    train_r2 = r2_score(train_true_list, train_pred_list)
del train_pool_x_t, train_pool_y_t, y_pred_train

# --- 4. Evaluate ONCE on the truly held-out test set (never touched above) ---
held_out_x_t = torch.tensor(held_out_x, dtype=torch.float32).unsqueeze(1).to(device)
held_out_y_t = torch.tensor(held_out_y, dtype=torch.float32).unsqueeze(1).to(device)
final_model.eval()
with torch.no_grad():
    y_pred_held = final_model(held_out_x_t)
    pred_list = [v[0] for v in y_pred_held.cpu().numpy().tolist()]
    true_list = [v[0] for v in held_out_y_t.cpu().numpy().tolist()]
    test_r2 = r2_score(true_list, pred_list)

del final_model, inner_train_x_t, inner_train_y_t, inner_val_x_t, inner_val_y_t, held_out_x_t, held_out_y_t
torch.cuda.empty_cache()
gc.collect()

File = open(args.output_file, 'a')
File.write(f"{TxID}\t{mean_r2}\t{mean_cv_r2}\t{train_r2}\t{test_r2}\t{len(train_pool_df)}\t{len(held_out_df)}\n")
File.flush()
