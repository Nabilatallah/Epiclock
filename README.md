# Epiclock
 Epigenetic clocks, based on DNA methylation patterns, is a machine learning based powerful predictor of biological age and health outcomes.
# 🧬 Epiclock: Machine Learning-Based Epigenetic Clock

## Overview
Epiclock is a complete, reproducible pipeline for building and interpreting a **DNA methylation-based epigenetic clock** using the publicly available **GSE40279** dataset.  
The pipeline includes:

- 📥 **Automated Data Retrieval** from GEO
- 🧹 **Preprocessing & Alignment** of methylation beta values and sample metadata
- 🤖 **Multiple Machine Learning Models** (Random Forest, ElasticNet, Ridge, Lasso, XGBoost, LightGBM, SVR)
- 📊 **Model Interpretation** using SHAP values
- 🔍 **Biological Annotation** of CpGs
- 📈 **Gene Set Enrichment & Visualization**

---

## Project Structure
epiclock/
├── code/ # All scripts
│   
├── data/ # Downloaded datasets
├── results/ # Processed outputs & logs
├── img/ # Plots and figures


---

## Installation

### 1️⃣ Clone the Repository
```bash
git clone https://github.com/Nabilatallah/epiclock.git
cd epiclock

If you use this pipeline in your research, please cite:

Atallah, N., Epiclock: A Machine Learning-Based Epigenetic Clock Pipeline Using DNA Methylation Data, 2025.

License
MIT License — see LICENSE file for details.

Contact
👤 Nabil Atallah
📧 n.atallah@northeastern.edu
🏛 Northeastern University, Boston, MA




