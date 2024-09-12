# Reconstructing the diversity dynamics of Late Quaternary grasslands through application of deep learning on superresolution images of fossil Poaceae pollen

<p align="center">
  <img src="https://github.com/paleopollen/Pollen_Diversity_Dynamics/blob/main/Figures_Rutundu/Figure_2_MountKenya_Map.png" width="950" title="hover text">
</p>

<p align="center">
  <img src="https://github.com/paleopollen/Pollen_Diversity_Dynamics/blob/main/Figures_Rutundu/Figures_1-4.png" width="750" title="hover text">
</p>
  
# Abstract
Despite its abundance in the fossil record, grass pollen is largely overlooked as a source of ecological and evolutionary data because species and genera cannot be easily discriminated visually. However, superresolution imaging and deep learning can identify morphological differences among grass taxa by focusing on small variations in grain morphology and surface ornamentation. Using a semi-supervised learning strategy, we trained convolutional neural networks (CNNs) on image data from extant Poaceae species and unlabeled fossil samples. Semi-supervised learning improved the CNN models' capability to generalize feature recognition in fossil pollen specimens. Our models successfully capture morphological diversity among our 60 modern species and between C<sub>3</sub> and C<sub>4</sub> grasses. We applied our trained models to fossil pollen images from a 25,000-year sediment core from Lake Rutundu, Mt Kenya, and identified a strong correlation between shifts in grass diversity and atmospheric CO<sub>2</sub> levels, temperature, precipitation, and fire frequency as measured by charcoal abundance. We quantified grass diversity by focusing on morphotypic variability, calculating Shannon entropy and morphotypic counts from the probability density function of specimens’ CNN features for each core depth. Predicted C<sub>3</sub> : C<sub>4</sub> abundances suggest a gradual increase in C<sub>4</sub> grass species correlated with rising temperature and fire activity. Our results demonstrate that machine-learned morphological features can significantly advance palynological analysis, enabling the estimation of biodiversity and distinction between C<sub>3</sub> and C<sub>4</sub> grass pollen using morphology alone.

# Significance Statement 
Deep learning and superresolution imaging are capable of solving some of the most intractable identification problems in fossil pollen analysis. The pollen of grass species are morphologically indistinguishable under traditional light and difficult to discriminate. However, superresolution imaging and deep learning successfully distinguishes the pollen of grass species. Features derived from convolutional neural networks quantify the biological and physiological diversity of grass pollen assemblages and can be applied without a priori knowledge of the species present, allowing reconstructions of changes in grass diversity and C<sub>4</sub> abundance. This approach unlocks new ecological information preserved in the abundant grass pollen record.

# Main Structure 
There are three main folders in this repository:
1. [Training and Classification](https://github.com/paleopollen/Pollen_Diversity_Dynamics/tree/main/00_Training_and_Classification): Scripts for training the two classification models described in the paper using two modalities: maximum intensity projection (MIP) images, and patches.
2. [Biodiversity Estimation](https://github.com/paleopollen/Pollen_Diversity_Dynamics/tree/main/01_Diversity_Estimation): Scripts for running the ecological simulations described in the paper and for applying Shannon entropy to calculate morphological diversity along the Lake Rutundu sediment core over the past 25,000 years. 
3. [Photosynthetic Pathway Analysis](https://github.com/paleopollen/Pollen_Diversity_Dynamics/tree/main/02_Photosynthetic_Pathway_Analysis): Scripts for detecting morphological differences between $C_3$ and $C_4$ grass pollen while accounting for phylogenetic relatedness, and for developing a random forest classifier to identify the photosynthetic pathway ($C_3$ or $C_4$) of grass fossil pollen based on morphology alone. 

# Hardware Specifications
Experiments were conducted on an NVIDIA GeForce RTX3090 GPU card with 24 GB of memory and an NVIDIA A100 SXM4 card with 40 GB of memory. We used the [PyTorch toolbox](https://pytorch.org/) for training neural networks. Additionally, some analyses were performed using R on a standard CPU. 

