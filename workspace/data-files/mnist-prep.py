import torch
from torchvision import datasets, transforms
import matplotlib.pyplot as plt
import numpy as np
import json
import random

def download_and_process_mnist(digit_to_extract):
    # Define a transform to normalize the data, then compress pixel values
    transform = transforms.Compose([
        transforms.ToTensor(),  # Converts to [0, 1] range
        transforms.Lambda(lambda x: torch.floor(x * 16))  # Compress to [0, 16]
    ])

    # Download MNIST data
    train_data = datasets.MNIST(root='./data', train=True, download=True, transform=transform)
    test_data = datasets.MNIST(root='./data', train=False, download=True, transform=transform)

    # Filter for the specific digit
    train_data_filtered = [(img, label) for img, label in train_data if label == digit_to_extract]
    test_data_filtered = [(img, label) for img, label in test_data if label == digit_to_extract]

    # Combine filtered datasets
    combined_data = train_data_filtered + test_data_filtered

    # Prepare JSON structure
    images_json = []

    for img, label in combined_data:
        # Flatten the image to a single list and convert to numpy array then to list for JSON serialization
        images_json.append(img.numpy().squeeze().astype(int).flatten().tolist())

    # Save to JSON file
    json_file_path = 'mnist-compressed.json'
    with open(json_file_path, 'w') as f:
        json.dump(images_json, f)
    
    print(f'MNIST data for digit {digit_to_extract} saved to {json_file_path}')

    # Visualize a random sample of digits
    fig, axes = plt.subplots(1, 10, figsize=(15, 1.5))
    random_samples = random.sample(combined_data, 10)
    for ax, (img, label) in zip(axes, random_samples):
        ax.imshow(img.squeeze(), cmap='gray', vmin=0, vmax=16)
        ax.set_title(label)
        ax.axis('off')

    plt.tight_layout()
    plt.show()

# Specify the digit you want to extract
digit_to_extract = 3  # For example, extracting images for the digit 2

# Call the function to process and visualize MNIST data for the specified digit
download_and_process_mnist(digit_to_extract)

