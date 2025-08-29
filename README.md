# net2rank

A new methodology that combines analysis of disease-specific omics data, network-based protein embeddings, and supervised machine learning to map the disease-protein associations.

## Citation


## Installation

```bash
git clone https://github.com/larsjuhljensen/net2rank.git

#install uv
pip install uv

#install net2rank
cd net2rank
uv sync
uv pip install .

# make sure you activate the virtual environment
source .venv/bin/activate

# double check
which python
```

## Cross validation
```bash
python main.py cross_validation \
--train_file data/train/aortic_aneurysm.olink.tsv \
--file_type label

python main.py cross_validation \
--train_file data/train/colorectal_adenocarcinoma.mutations.intogen.tsv \
--file_type list

python main.py cross_validation \
--train_file data/train/melanoma.mutations.intogen.tsv \
--file_type list

python main.py cross_validation \
--train_file data/train/diffuse_large_b-cell_lymphoma.mutations.intogen.tsv \
--file_type list

python main.py cross_validation \
--train_file data/train/atopic_dermatitis.integrated.tsv \
--file_type pvalue \
--pos_size 1000

python main.py cross_validation \
--train_file data/train/ulcerative_colitis.integrated.tsv \
--file_type pvalue \
--pos_size 1300

python main.py cross_validation \
--train_file data/train/focal_epilepsy.rnaseq.kjaer_guelfi_consensus.tsv \
--file_type list
```

## Train and test
```bash
python main.py train_test \
--train_file data/train/atopic_dermatitis.integrated.tsv \
--file_type pvalue \
--test_file data/test/atopic_dermatitis.gold_standard.balanced.tsv \
--pos_size 1000

python main.py train_test \
--train_file data/train/ulcerative_colitis.integrated.tsv \
--file_type pvalue \
--test_file data/test/ulcerative_colitis.gold_standard.balanced.tsv \
--pos_size 1300

python main.py train_test \
--train_file data/train/focal_epilepsy.rnaseq.kjaer_guelfi_consensus.tsv \
--file_type list \
--test_file data/test/focal_epilepsy.gold_standard.balanced.tsv

python main.py train_test \
--train_file data/train/colorectal_adenocarcinoma.mutations.intogen.tsv \
--file_type list \
--test_file data/test/colorectal_adenocarcinoma.gold_standard.balanced.tsv

python main.py train_test \
--train_file data/train/melanoma.mutations.intogen.tsv \
--file_type list \
--test_file data/test/melanoma.gold_standard.balanced.tsv

python main.py train_test \
--train_file data/train/diffuse_large_b-cell_lymphoma.mutations.intogen.tsv \
--file_type list \
--test_file data/test/diffuse_large_b-cell_lymphoma.gold_standard.balanced.tsv

python main.py train_test \
--train_file data/train/aortic_aneurysm.olink.tsv \
--file_type label \
--test_file data/test/aortic_aneurysm.gold_standard.balanced.tsv

```

## Network visualization
We used Cytoscape stringAPP to visualize newtworks, and used py4cytoscape for some automation. Please check the notebook `notebooks/enrichment.ipynb` for details. We also provided our Cytoscape session here: https://zenodo.org/records/16919169