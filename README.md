# Deep learning of fossil pollen morphology reveals 25,000 years of ecological change in East African grasslands

<p align="center">
  <img src="https://github.com/paleopollen/Pollen_Diversity_Dynamics/blob/main/Figures_Rutundu/Rutundu_Figure_2_Revised.png" width="800" title="hover text">
</p>

<br><br>

<p align="center">
  <img src="https://github.com/paleopollen/Pollen_Diversity_Dynamics/blob/main/Figures_Rutundu/Rutundu_Figure_1_4_Revised.png" width="750" title="hover text">
</p>

<br><br>

<p align="center">
  <img src="https://github.com/paleopollen/Pollen_Diversity_Dynamics/blob/main/Figures_Rutundu/Rutundu_Figure_5_Revised.png" width="750" title="hover text">
</p>

# Abstract
Grass pollen is largely overlooked in investigations of grassland evolution because the pollen of most species cannot be differentiated using traditional optical microscopes. However, superresolution microscopy captures distinct morphological variations across grass pollen species and deep learning can quantify such variations. Using a semi-supervised deep learning strategy, we trained convolutional neural networks (CNNs) on superresolution pollen images of modern grasses and unlabeled fossil Poaceae. The learned features reflected both the taxonomic diversity of grass communities along an elevational gradient and morphological differences between C<sub>3</sub> and C<sub>4</sub> species. We applied our trained models to fossil grass pollen assemblages from a 25,000-year lake-sediment record from eastern equatorial Africa (Mt. Kenya) and correlated past shifts in grass diversity with atmospheric CO<sub>2</sub> concentration and proxy records of local temperature, precipitation, and fire occurrence. We quantified changes in grass diversity using morphological variability of fossil pollen assemblages, approximated by the Shannon entropy of CNN features. Our data show that grassland species diversity decreased substantially between 21,500 and 16,000 years ago, coincident with most severe regional cooling during the last ice age. C<sub>3</sub>:C<sub>4</sub> ratios reconstructed using a gradient boosted decision tree classifier infer a gradual decrease in C<sub>4</sub> grasses since the late-glacial to Holocene transition, associated with decreasing fire activity and elevated temperatures. Our results demonstrate that CNN features of pollen morphology can advance palynological analysis, enabling robust estimation of grass diversity and C<sub>3</sub>:C<sub>4</sub> ratios in ancient grassland ecosystems.

# Significance Statement 
Grasses are an ecologically important and diverse plant clade whose pollen, being morphologically indistinguishable under standard optical microscopy, has historically provided limited insight into the ecological history of grasslands. We show, however, that grass pollen can be differentiated through deep learning analyses of superresolution images. Abstracted features derived from deep neural networks can be used to quantify the biological and physiological diversity of grass pollen assemblages, without a priori knowledge of the species present, and reconstruct past changes in grass diversity and the relative abundance of C<sub>4</sub> grasses in ancient grasslands. This approach unlocks previously unattainable ecological information from fossil pollen and demonstrates that deep learning can solve seemingly intractable identification problems and advance the reconstruction of past vegetation dynamics. 

# Main Structure 
There are three main folders in this repository:
1. [Training and Classification](https://github.com/paleopollen/Pollen_Diversity_Dynamics/tree/main/00_Training_and_Classification): Scripts for training the two classification models described in the paper using two modalities: maximum intensity projection (MIP) images, and patches.
2. [Species Diversity Estimation](https://github.com/paleopollen/Pollen_Diversity_Dynamics/tree/main/01_Species_Diversity_Estimation): Scripts for running the ecological simulations described in the paper and for estimating morphological diversity in pollen assemblages using KDE entropy on CNN-derived features, without requiring species identification. Methods can be applied to reconstruct diversity changes along sediment core (e.g., Lake Rutundu over the past ~25,000 years). 
3. [Photosynthetic Pathway Identification](https://github.com/paleopollen/Pollen_Diversity_Dynamics/tree/main/02_Photosynthetic_Pathway_Identification): Scripts for detecting morphological differences between C<sub>3</sub> and C<sub>4</sub> grass pollen while accounting for phylogenetic relatedness, and for developing a random forest classifier to identify the photosynthetic pathway (C<sub>3</sub> or C<sub>4</sub>) of grass fossil pollen based on morphology alone. 

# Hardware Specifications
Experiments were conducted on an NVIDIA GeForce RTX3090 GPU card with 24 GB of memory and an NVIDIA A100 SXM4 card with 40 GB of memory. We used the [PyTorch toolbox](https://pytorch.org/) for training neural networks. Additionally, some analyses were performed using R on a standard CPU.
