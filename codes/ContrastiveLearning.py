#### Contrastive Learning for determining significant features in CCCs ####

import sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.linear_model import LogisticRegression


PositiveData = pd.read_csv(sys.argv[1], sep=',', header=None)
NegativeData = pd.read_csv(sys.argv[2], sep=',', header=None)
FullData = pd.read_csv(sys.argv[3], sep=',', header=None, index_col=0)
row_names = FullData.index.astype(str)
X = np.vstack([PositiveData.values, NegativeData.values])
y = np.array([1]*len(PositiveData) + [0]*len(NegativeData))
data = torch.tensor(X, dtype=torch.float32)
fulldata = torch.tensor(FullData.values, dtype=torch.float32)
labels = torch.tensor(y, dtype=torch.float32)

if torch.cuda.is_available():
    device = torch.device('cuda')
    print('Using GPU:', torch.cuda.get_device_name(0))
else:
    device = torch.device('cpu')
    print('Using CPU')

data = data.to(device)
labels = labels.to(device)

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
    labels = labels.view(-1, 1)
    mask = torch.eq(labels, labels.T).float().to(device)
    features = F.normalize(features, dim=1)
    anchor_dot_contrast = torch.div(torch.matmul(features, features.T), temperature)
    logits_max, _ = torch.max(anchor_dot_contrast, dim=1, keepdim=True)
    logits = anchor_dot_contrast - logits_max.detach()
    logits_mask = torch.ones_like(mask) - torch.eye(mask.shape[0], device=device)
    mask = mask * logits_mask
    exp_logits = torch.exp(logits) * logits_mask
    log_prob = logits - torch.log(exp_logits.sum(1, keepdim=True) + 1e-12)
    mean_log_prob_pos = (mask * log_prob).sum(1) / (mask.sum(1) + 1e-12)
    loss = -mean_log_prob_pos.mean()
    return loss

input_dim = data.shape[1]
encoder = Encoder(input_dim).to(device)
optimizer = torch.optim.Adam(encoder.parameters(), lr=1e-3)
epochs = 300
batch_size = 64

for epoch in range(epochs):
    perm = torch.randperm(data.size(0))
    data_shuffled = data[perm]
    labels_shuffled = labels[perm]
    for i in range(0, data.size(0), batch_size):
        x_batch = data_shuffled[i:i+batch_size]
        y_batch = labels_shuffled[i:i+batch_size]
        if x_batch.size(0) < 2:
            continue
        z = encoder(x_batch)
        loss = supcon_loss(z, y_batch)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    if (epoch+1) % 10 == 0:
        print(f"Epoch {epoch+1}, SupCon Loss: {loss.item():.4f}")


train_features = encoder(data).detach().cpu().numpy()
train_labels = labels.cpu().numpy()
clf = LogisticRegression(max_iter=1000)
clf.fit(train_features, train_labels)

fulldata_features = encoder(fulldata.to(device)).detach().cpu().numpy()
fulldata_pred = clf.predict(fulldata_features)

pd.DataFrame({'CCC_pair': row_names, 'prediction': fulldata_pred}).to_csv('fulldata_prediction.csv', index=False)
