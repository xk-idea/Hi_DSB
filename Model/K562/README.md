# Overview
This project leverages Hi-C data to construct genomic graphs, integrates multi-omics data, and employs graph contrastive learning models to predict Double-Strand Break (DSB) sites across the entire genome. By utilizing advanced graph-based methodologies, this pipeline offers a comprehensive approach to understanding genomic instability and its implications in various biological processes.

# Installation
You can use the following command to install the necessary environment:
  `pip install -r requirements.txt`

# GraphMaker
To build the primary genomic graph using Hi-C data:
  `jupyter notebook 1.Graph_maker.ipynb`
This notebook processes Hi-C datasets to generate the graph structure essential for downstream analysis.

# Feature-Similarity Graph (Optional)
If you wish to incorporate feature similarity into the graph construction:
  `jupyter notebook 2.Sim2.ipynb`

# Model training
The operation instructions are:  
  `python GraphstJK_Cross-validation_4feats.py`  
The training script includes built-in cross-validation and parameter settings to facilitate immediate usage.

# Model explain
Interpret the trained model to understand edge contributions: 
  `python EX_edge.py`  

# Results
Upon completion of model training and interpretation, all trained models and execution results are automatically saved in the `results` directory. This includes:
  1.Trained model checkpoints
  2.Model prediction performance evaluation
  3.Interpretation reports and visualizations
