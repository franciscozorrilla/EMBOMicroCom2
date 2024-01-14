# 🦠 Systems Biology Practical: FBA with GEMs 🧬

[![DOI](https://zenodo.org/badge/740951705.svg)](https://zenodo.org/doi/10.5281/zenodo.10475949)

## 📜 About

This repo hosts training materials for the University of Cambridge [Part III Systems Biology course](https://www.sysbiol.cam.ac.uk/Part%20III). More specifically, these are practical exercises for the computational session on flux balance analysis (FBA) using genome-scale metabolic models (GEMs). Some of the text and code chunks in the following exercises are based on selected portions of the [cobrapy documentation](https://cobrapy.readthedocs.io/en/latest/).

### Learning outcomes

#### [Exercise 1](https://github.com/franciscozorrilla/systems-biology-fba-practical/blob/main/notebooks/1_fba.ipynb)
- **1.1**: Import a metabolic reconstruction
- **1.2**: Inspect the reactions of your model
- **1.3**: Inspect the metabolites in your model
- **1.4**: Inspect the genes in your model
- **1.4.1**: Perform in-silico gene knockout experiments

#### [Exercise 2](https://github.com/franciscozorrilla/systems-biology-fba-practical/blob/main/notebooks/2_fba.ipynb)
- **2.1**: Modify growth medium of your reconstruction
- **2.2**: Perform gene essentiality analysis under different conditions
- **2.3**: Case study, simulate the Crabtree effect in yeast

## 🚚 Software requirements

The following packages will already be installed on the provided machines:

* [cobrapy](https://opencobra.github.io/cobrapy/)
* [matplotlib](https://matplotlib.org/stable/)
* [numpy](https://numpy.org/install/)
* [notebook](https://jupyter.org/install#jupyter-notebook)

To install on your own computers use the [pip package manager](https://pip.pypa.io/en/stable/getting-started/).

Using the `requirements.txt` files:

```
$ pip install -r requirements.txt
```

Or from scratch:

```
$ pip install cobra matplotlib numpy notebook
```

## 🛠️ Usage

You will find everything you need for these practical exercises under the `/notebooks` folder. Each tutorial is uploaded as `.html`, `.md`, `.tex`, and `.ipynb` files. As you carry out the exercises using the the python notebooks, you will find the original text and results in the html, markdown, and latex files. 

### Running locally

You may run these exercises on your local machine by cloning this repo:

```
$ git clone https://github.com/franciscozorrilla/systems-biology-fba-practical.git
```

Navigate into cloned repo folder:

```
$ cd systems-biology-fba-practical
```

Launch interactive jupyter notebook session:

```
$ jupyter notebook
```

This should launch a browser window, where you can click on the `/notebooks` folder and start the exercises.

### Running on the cloud 

Alternatively launch an interactive binder session by clicking below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/franciscozorrilla/systems-biology-fba-practical/HEAD)

## 🧠 Further exercises

* [metaGEM](https://github.com/franciscozorrilla/unseenbio_metaGEM): Metagenomics-driven metabolic modelling tutorial using unseen bio samples
* [SymbNET](https://github.com/franciscozorrilla/SymbNET): From metagenomics to metabolic Interactions 
* [EMBOMicroCom](https://github.com/franciscozorrilla/EMBOMicroCom): Metabolite and species dynamics in microbial communities

## 👷 Contributors

* Originally developed by Arianna Basile and Kiran Patil in 2023
* Updated by Francisco Zorrilla and Arianna Basile in 2024


## ☎️ Contact

Feel free to get in touch with us by raising an issue or sending an e-mail to Arianna Basile (ab2851@mrc-tox.cam.ac.uk) or Francisco Zorrilla (fz274@cam.ac.uk).
