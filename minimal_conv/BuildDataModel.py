import torch
from torch.utils.data import Dataset,DataLoader
import os   
#from torchvision import transforms
import pandas as pd
import nibabel as nib
import numpy as np
import math
######################################################################## Custom dataset class##########################################################################################
class CustomDataset(Dataset):
    def __init__(self, csv_file, root_dir, transform=None):
        self.annotations = pd.read_csv(csv_file)
        self.root_dir = root_dir
        self.transform = transform
    def __len__(self):
        return len(self.annotations)
    def __getitem__(self, index):
        img_path = os.path.join(self.root_dir, self.annotations.iloc[index, 0])
        image_fname = self.annotations.iloc[index, 0]
        image = nib.load(img_path)
        image = image.get_fdata()
        y_label = torch.tensor(int(self.annotations.iloc[index, 1]))
        if self.transform:
            image = self.transform(image)
        return (image, y_label, image_fname)

# Custom transform class to read in nifti file, apply mask, then convert to tensor
class ToTensor:
    def __call__(self, image):
        # Apply mask
        mask = nib.load('/home/users/hyang336/retro_timing_fMRI/ROIs/Schaefer 2018 400 Parcels 17 Networks/MNI152NLin2009cAsym_res-02_atlas-Schaefer2018_desc-400Parcels17Networks_lSomMotA_bVisCent_ExStr.nii')
        mask = mask.get_fdata()
        image = np.nan_to_num(image)#convert nan to 0
        image = np.multiply(image, mask)
        # Convert to tensor
        image = torch.from_numpy(image)
        image = image.unsqueeze(0)#add channel dimension, which is 1 in this case
        return image

###################################################################### Data loader######################################################################################################
# Set device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Hyperparameters
input_size = 97*115*97
learning_rate = 0.001
batch_size = 1
num_epochs = 10

# Load Data
dataset = CustomDataset(csv_file='/scratch/users/hyang336/retro_timing_fMRI_data/yaam_torch_data/annotations.csv', root_dir='/scratch/users/hyang336/retro_timing_fMRI_data/yaam_torch_data/', transform=ToTensor())
datalength=len(dataset)
train_set, test_set = torch.utils.data.random_split(dataset, [math.ceil(0.7*datalength), math.floor(0.3*datalength)])
train_loader = DataLoader(dataset=train_set, batch_size=batch_size, shuffle=True)
test_loader = DataLoader(dataset=test_set, batch_size=batch_size, shuffle=True)

# Model
class mininal3DConvNet(torch.nn.Module):
    def __init__(self):
        super(mininal3DConvNet, self).__init__()
        self.conv_stack=torch.nn.Sequential(
            torch.nn.Conv3d(1, 8, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.MaxPool3d(kernel_size=2, stride=2),
            torch.nn.Conv3d(8, 32, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.MaxPool3d(kernel_size=2, stride=2),
            torch.nn.Conv3d(32, 128, kernel_size=3, stride=1, padding=1),
            torch.nn.ReLU(),
            torch.nn.MaxPool3d(kernel_size=2, stride=2)
        )
        self.fc_stack=torch.nn.Sequential(
            torch.nn.Linear(128*12*14*12, 512),
            torch.nn.ReLU(),
            torch.nn.Linear(512, 1)
        )

    def forward(self, x):
        x = self.conv_stack(x)
        x = torch.flatten(x,1)
        x = self.fc_stack(x)
        return x

# Create model
model=mininal3DConvNet().to(device)

# Loss and optimizer
loss_fn=torch.nn.MSELoss()
optimizer=torch.optim.Adam(model.parameters(), lr=learning_rate)

# Train and test functions
def train(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    for batch, (X, y, z) in enumerate(dataloader):
        X, y = X.to(device,dtype=torch.float), y.to(device,dtype=torch.float)

        # Compute prediction error
        pred = model(X)
        y = y.view_as(pred)  # Reshape y to have the same size as pred
        loss = loss_fn(pred, y)

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if batch % 100 == 0:
            loss, current = loss.item(), batch * len(X)
            print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
            print(f"image_fname: {z}  [{current:>5d}/{size:>5d}]")

def test(dataloader, model, loss_fn):
    model.eval()
    size = len(dataloader.dataset)
    test_loss, correct = 0, 0

    with torch.no_grad():
        for X, y, z in dataloader:
            X, y = X.to(device,dtype=torch.float), y.to(device,dtype=torch.float)
            pred = model(X)
            y = y.view_as(pred)  # Reshape y to have the same size as pred
            test_loss += loss_fn(pred, y).item()
            correct += (pred.argmax(1) == y).type(torch.float).sum().item()

    test_loss /= size
    correct /= size
    print(f"Test Error: \n Accuracy: {(100*correct):>0.1f}%, Avg loss: {test_loss:>8f} \n")


for t in range(num_epochs):
    print(f"Epoch {t+1}\n-------------------------------")
    train(train_loader, model, loss_fn, optimizer)
    test(test_loader, model, loss_fn)

torch.save(model.state_dict(), "/scratch/users/hyang336/retro_timing_fMRI_data/minimalConv/model.pth")    
print("Done!")
