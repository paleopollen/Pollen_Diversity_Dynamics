# Biodiversity dynamics of Late Quaternary East African grasslands reconstructed using fossil Poaceae pollen and deep learning 

<p align="center">
  <img src="https://github.com/madaime2/Pollen_Biodiversity_Reconstruction/blob/main/Grass_Project_Figures/Fig1_Grass_Pollen_Project.png" width="950" title="hover text">
</p>

<div style="display: flex; justify-content: center; align-items: center;">
  <img src="https://github.com/madaime2/Pollen_Biodiversity_Reconstruction/blob/main/Grass_Project_Figures/Fig2_Grass_Pollen_Project.png" width="350" title="hover text">
  
  <img src="https://github.com/madaime2/Pollen_Biodiversity_Reconstruction/blob/main/Grass_Project_Figures/Fig3_Grass_Pollen_Project.png" width="350" title="hover text" style="margin: 0 15px;">

  <img src="https://github.com/madaime2/Pollen_Biodiversity_Reconstruction/blob/main/Grass_Project_Figures/Fig4_Grass_Pollen_Project.png" width="350" title="hover text">
</div>
  

# Abstract
Despite its abundance in the fossil record, grass pollen is largely overlooked as a source of ecological and evolutionary data because species and genera cannot be easily discriminated visually. However, superresolution imaging and deep learning can identify morphological differences among grass taxa by focusing on small variations in grain morphology and surface ornamentation. Using a semi-supervised learning strategy, we trained convolutional neural networks (CNNs) on image data from extant Poaceae species and unlabeled fossil samples. Semi-supervised learning improved the CNN models' capability to generalize feature recognition in fossil pollen specimens. Our models successfully capture morphological diversity among our 60 modern species and between $C_3$ and $C_4$ grasses. We applied our trained models to fossil pollen images from a 25,000-year sediment core from Lake Rutundu, Mt Kenya, and identified a strong correlation between shifts in grass diversity and atmospheric $CO_2$ levels, temperature, precipitation, and fire frequency as measured by charcoal abundance. We quantified grass diversity by focusing on morphotypic variability, calculating Shannon entropy and morphotypic counts from the probability density function of specimens’ CNN features for each core depth. Predicted $C_3$ : $C_4$ abundances suggest a gradual increase in $C_4$ grass species correlated with rising temperature and fire activity. Our results demonstrate that machine-learned morphological features can significantly advance palynological analysis, enabling the estimation of biodiversity and distinction between $C_3$ and $C_4$ grass pollen using morphology alone.

# Significance Statement 
Deep learning and superresolution imaging are capable of solving some of the most intractable identification problems in fossil pollen analysis. The pollen of grass species are morphologically indistinguishable under traditional light and difficult to discriminate. However, superresolution imaging and deep learning successfully distinguishes the pollen of grass species. Features derived from convolutional neural networks quantify the biological and physiological diversity of grass pollen assemblages and can be applied without a priori knowledge of the species present, allowing reconstructions of changes in grass diversity and $C_4$ abundance. This approach unlocks new ecological information preserved in the abundant grass pollen record.

# Main Structure 
There are four folders in this repository:

# Hardware Specifications
Experiments were conducted on an NVIDIA GeForce RTX3090 GPU card with 24 GB of memory and an NVIDIA A100 SXM4 card with 40 GB of memory. We used the [PyTorch toolbox](https://pytorch.org/) for training neural networks.

