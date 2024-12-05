import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torchvision.models import (resnet18, ResNet18_Weights,
                                resnet34, ResNet34_Weights,
                                resnet50, ResNet50_Weights,
                                resnet101, ResNet101_Weights, 
                                resnet152, ResNet152_Weights,
                                resnext50_32x4d,
                                resnext101_32x8d,
                                resnext101_64x4d)

class C2P1_FC1(nn.Module):
    """
    Two convolutional layers followed by one pooling layer and one fully connected layer.

    Args
    ----
    nchannels : int
        Number of energy bins.
    nphases : int
        Number of phase bins.
    out_channels : int, optional
        Number of channels in the first convolutional layer. Default is 16.
    kernel_size : int, optional
        Kernel size of the convolutional layers. Default is 3.
    pool_kernel_size : int, optional
        Kernel size of the pooling layer. Default is 2.
    pool_stride : int, optional
        Stride of the pooling layer. Default is 2.
    fc_out_features : int, optional
        Number of features in the output of the fully connected layer. Default is 10.

    Attributes
    ----------
    conv1 : torch.nn.LazyConv2d
        First convolutional layer.
    conv2 : torch.nn.LazyConv2d
        Second convolutional layer.
    pool : torch.nn.MaxPool2d
        Pooling layer.
    fc : torch.nn.LazyLinear
        Fully connected layer.

    """
    def __init__(self, nchannels, nphases,
                 out_channels=16, kernel_size=3, 
                 pool_kernel_size=2, pool_stride=2, 
                 fc_out_features=10):
        super().__init__()

        self.conv1 = nn.LazyConv2d(out_channels=out_channels, kernel_size=kernel_size, padding=kernel_size//2)
        self.conv2 = nn.LazyConv2d(out_channels=out_channels*2, kernel_size=kernel_size, padding=kernel_size//2)
        self.pool = nn.MaxPool2d(kernel_size=pool_kernel_size, stride=pool_stride)
        self.fc = nn.LazyLinear(out_features=fc_out_features)
        self.nchannels = nchannels
        self.nphases = nphases        

    def forward(self, x):
        x = x.view(-1, 1, self.nchannels, self.nphases)
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = self.pool(x)
        x = torch.flatten(x, 1)
        x = F.relu(self.fc(x))

        return x


class C2P1_C1P1_FC1(nn.Module):
    """
    Two convolutional layers followed by one pooling layer, followed by one convolutional layer, 
    followed by one pooling layer, and one fully connected layer.

    Args
    ----
    nchannels : int
        Number of energy bins.
    nphases : int
        Number of phase bins.
    out_channels : int, optional
        Number of channels in the first convolutional layer. Default is 16.
    kernel_size : int, optional
        Kernel size of the convolutional layers. Default is 3.
    pool_kernel_size : int, optional
        Kernel size of the pooling layers. Default is 2.
    pool_stride : int, optional
        Stride of the pooling layers. Default is 2.
    fc_out_features : int, optional
        Number of features in the output of the fully connected layer. Default is 10.

    Attributes
    ----------
    conv1 : torch.nn.LazyConv2d
        First convolutional layer.
    conv2 : torch.nn.LazyConv2d
        Second convolutional layer.
    pool1 : torch.nn.MaxPool2d
        First pooling layer.
    conv3 : torch.nn.LazyConv2d
        Third convolutional layer.
    pool2 : torch.nn.MaxPool2d
        Second pooling layer.
    fc : torch.nn.LazyLinear
        Fully connected layer.

    """
    def __init__(self, nchannels, nphases,
                 out_channels=16, kernel_size=3, 
                 pool_kernel_size=2, pool_stride=2, 
                 fc_out_features=10):
        super().__init__()

        self.conv1 = nn.LazyConv2d(out_channels=out_channels, kernel_size=kernel_size, padding=kernel_size//2)
        self.conv2 = nn.LazyConv2d(out_channels=out_channels*2, kernel_size=kernel_size, padding=kernel_size//2)
        self.pool1 = nn.MaxPool2d(kernel_size=pool_kernel_size, stride=pool_stride)
        self.conv3 = nn.LazyConv2d(out_channels=out_channels*4, kernel_size=kernel_size, padding=kernel_size//2)
        self.pool2 = nn.MaxPool2d(kernel_size=pool_kernel_size, stride=pool_stride)
        self.fc = nn.LazyLinear(out_features=fc_out_features)
        self.nchannels = nchannels
        self.nphases = nphases    

    def forward(self, x):
        x = x.view(-1, 1, self.nchannels, self.nphases)
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = self.pool1(x)
        x = F.relu(self.conv3(x))
        x = self.pool2(x)
        x = torch.flatten(x, 1)
        x = F.relu(self.fc(x))

        return x
    
class embedding_resnet(nn.Module):
    """
    A class that wraps a pre-trained `ResNet/ResNeXt model <https://pytorch.org/vision/stable/models.html#general-information-on-pre-trained-weights>` 
    and unfreezes the last backbone and linear layers for fine-tuning.

    Args
    ----
    resnet_version (str): The version of ResNet/ResNeXt to use
        (one of 'resnet18', 'resnet34', 'resnet50',
        'resnet101', 'resnet152', 'resnext50',
        'resnext101', 'resnext101_64x4d').

    Attributes
    ----------
    resnet: The pre-trained ResNet/ResNeXt model
    fc: The linear layer at the end of the model
    """

    def __init__(self, resnet_version='resnet18'):
        super().__init__()
        if resnet_version == 'resnet18':
            self.resnet = resnet18(weights=ResNet18_Weights.DEFAULT)
            print('Using resnet18')
        elif resnet_version == 'resnet34':
            self.resnet = resnet34(weights=ResNet34_Weights.DEFAULT)
            print('Using resnet34')
        elif resnet_version == 'resnet50':
            self.resnet = resnet50(weights=ResNet50_Weights.DEFAULT)
            print('Using resnet50')
        elif resnet_version == 'resnet101':
            self.resnet = resnet101(weights=ResNet101_Weights.DEFAULT)
            print('Using resnet101')
        elif resnet_version == 'resnet152':
            self.resnet = resnet152(weights=ResNet152_Weights.DEFAULT)
            print('Using resnet152')
        elif resnet_version == 'resnext50':
            self.resnet = resnext50_32x4d(weights='DEFAULT')
        elif resnet_version == 'resnext101':
            self.resnet = resnext101_32x8d(weights='DEFAULT')
        elif resnet_version == 'resnext101_64x4d':
            self.resnet = resnext101_64x4d(weights='DEFAULT')
        else:
            raise ValueError(f'Invalid resnet version: {resnet_version}')

        # Freeze the backbone layers
        for param in self.resnet.parameters():
            param.requires_grad = False

        # Unfreeze the last few layers for fine-tuning
        for param in self.resnet.layer4.parameters():
            param.requires_grad = True
        for param in self.resnet.fc.parameters():
            param.requires_grad = True

        self.resnet.fc = nn.Linear(self.resnet.fc.in_features, 100)

    def forward(self, x):
        # Convert single-channel input to 3-channel
        x = x.unsqueeze(1)
        x = torch.cat([x, x, x], dim=1)
        x = self.resnet(x)
        return x
