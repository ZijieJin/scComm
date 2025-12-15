#### Contrastive Learning for determining significant features in CCCs ####
import sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans


if torch.cuda.is_available():
    device = torch.device('cuda')
    print('Using GPU:', torch.cuda.get_device_name(0))
elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
    device = torch.device('mps')
    print('Using Apple Silicon GPU')
else:
    device = torch.device('cpu')
    print('Using CPU')


class Encoder(nn.Module):
    def __init__(self, input_dim, hidden_dim=128, out_dim=64):
        super().__init__()
        self.fc1 = nn.Linear(input_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, out_dim)
    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return F.normalize(x, dim=1)


# SupCon Loss
def supcon_loss(features, labels, temperature=0.5):
    device = features.device
    labels = labels.view(-1, 1).long()
    mask = torch.eq(labels, labels.T).float().to(device)
    mask.fill_diagonal_(0)
    features = F.normalize(features, dim=1)
    anchor_dot_contrast = torch.div(torch.matmul(features, features.T), temperature)
    logits_max, _ = torch.max(anchor_dot_contrast, dim=1, keepdim=True)
    logits = anchor_dot_contrast - logits_max.detach()
    logits_mask = 1 - torch.eye(mask.shape[0], device=device)
    exp_logits = torch.exp(logits) * logits_mask
    log_prob = logits - torch.log(exp_logits.sum(1, keepdim=True) + 1e-12)
    denom = mask.sum(1)
    valid = denom > 0
    if valid.sum() == 0:
        return torch.tensor(0.0, device=device, requires_grad=True)
    mean_log_prob_pos = (mask * log_prob).sum(1) / (denom + 1e-12)
    loss = -mean_log_prob_pos[valid].mean()
    return loss



# Load Data
PositiveData = pd.read_csv('/Users/zijie/Desktop/CellCommunication/data/Revision1/R2.1.1/SimuData_2_8/Positive_3.csv', sep=',', index_col=0)
NegativeData = pd.read_csv('/Users/zijie/Desktop/CellCommunication/data/Revision1/R2.1.1/SimuData_2_8/Negative_3.csv', sep=',', index_col=0)
FullData = pd.read_csv('/Users/zijie/Desktop/CellCommunication/data/Revision1/R2.1.1/SimuData_2_8/FullData_3.csv', sep=',', index_col=0)

# PositiveData = pd.read_csv(sys.argv[1], sep=',', index_col=0)
# NegativeData = pd.read_csv(sys.argv[2], sep=',', index_col=0)
# FullData = pd.read_csv(sys.argv[3], sep=',', index_col=0)


# Prepare Data
common_cols = list(set(PositiveData.columns) & set(NegativeData.columns) & set(FullData.columns))
PositiveData = PositiveData[common_cols]
NegativeData = NegativeData[common_cols]
FullData = FullData[common_cols]
row_names = FullData.index.astype(str)
X = np.vstack([PositiveData.values, NegativeData.values])
y = np.array([1]*len(PositiveData) + [0]*len(NegativeData))
data = torch.tensor(X, dtype=torch.float32)
fulldata = torch.tensor(FullData.values, dtype=torch.float32)
labels = torch.tensor(y, dtype=torch.long)
data = data.to(device)
labels = labels.to(device)
input_dim = data.shape[1]
encoder = Encoder(input_dim).to(device)
optimizer = torch.optim.Adam(encoder.parameters(), lr=1e-4)
data_rs = data.sum(1)
mean_rs_pos = data_rs[labels==1].mean().item()
mean_rs_neg = data_rs[labels==0].mean().item()
labels[(data_rs < (mean_rs_pos + mean_rs_neg) / 2) & (labels == 1)] = 2
mean_rs_mid = data_rs[labels==2].mean().item()
labels[(data_rs < (mean_rs_pos + mean_rs_mid) / 2) & (labels == 1)] = 2

# Training Loop
encoder.train()
lastloss = 9999.0
epoch = 0
badcount = 0
while True:
    z = encoder(data)
    loss = supcon_loss(z, labels)
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    if (epoch+1) % 10 == 0:
        print(f"Epoch {epoch+1}, Avg SupCon Loss: {loss:.4f}")
        if loss.item() > lastloss - 1e-4 and epoch > 100:
            badcount += 1
            if badcount >= 5:
                print("Loss increased, stopping training.")
                break
        else:
            badcount = 0
            if epoch == 0 or loss.item() < lastloss:
                best_model_state = encoder.state_dict()
        lastloss = min(loss.item(), lastloss)
    epoch += 1
encoder.load_state_dict(best_model_state)

# Evaluation with Logistic Regression
encoder.eval()
with torch.no_grad():
    train_features = encoder(data).cpu().numpy()
train_labels = labels.cpu().numpy().astype(int)
clf = LogisticRegression(max_iter=1000)
clf.fit(train_features, train_labels)

# Predict on Full Data
with torch.no_grad():
    fulldata_features = encoder(fulldata.to(device)).cpu().numpy()
fulldata_pred = clf.predict(fulldata_features)

fulldata.sum(1)[fulldata_pred == 1]

# Save Predictions
fulldata_pred[fulldata_pred == 2] = 0
pd.DataFrame({'CCC_pair': row_names, 'prediction': fulldata_pred}).to_csv('fulldata_prediction.csv', index=False)
