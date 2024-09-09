# 🔬 [EMBOMicroCom2](https://www.embl.org/about/info/course-and-conference-office/events/mcd24-01/): Metabolite and species dynamics in microbial communities 🧬

# 💻 Metabolic modeling tutorial (Day 5)

## 📜 About

This repo hosts training materials for the second installment of EMBO's Practical Course [Metabolite and species dynamics in microbial communities](https://www.embl.org/about/info/course-and-conference-office/events/mcd24-01/). More specifically, these are practical exercises for the computational session on flux balance analysis (FBA) using genome-scale metabolic models (GEMs). Some of the text and code chunks in the following exercises are based on selected portions of the [cobrapy documentation](https://cobrapy.readthedocs.io/en/latest/).

### Learning outcomes

#### [Exercise 1](https://github.com/franciscozorrilla/EMBOMicroCom2/blob/main/notebooks/1_fba.ipynb)
- **1.1**: Import a metabolic reconstruction
- **1.2**: Inspect the reactions of your model
- **1.3**: Inspect the metabolites in your model
- **1.4**: Inspect the genes in your model
- **1.4.1**: Perform in-silico gene knockout experiments

#### [Exercise 2](https://github.com/franciscozorrilla/EMBOMicroCom2/blob/main/notebooks/2_fba.ipynb)
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
$ git clone https://github.com/franciscozorrilla/EMBOMicroCom2.git
```

Navigate into cloned repo folder:

```
$ cd EMBOMicroCom2
```

Launch interactive jupyter notebook session:

```
$ jupyter notebook
```

This should launch a browser window, where you can click on the `/notebooks` folder and start the exercises.

### Running on the cloud 

Alternatively launch an interactive binder session by clicking below:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/franciscozorrilla/EMBOMicroCom2/HEAD)

## 🧠 Further exercises

* [metaGEM](https://github.com/franciscozorrilla/unseenbio_metaGEM): Metagenomics-driven metabolic modelling tutorial using unseen bio samples
* [SymbNET](https://github.com/franciscozorrilla/SymbNET): From metagenomics to metabolic interactions 
* [EMBOMicroCom](https://github.com/franciscozorrilla/EMBOMicroCom): Metabolite and species dynamics in microbial communities
* [OLISSIPO](https://github.com/arianccbasile/FBA_OLISSIPO_Winter_School): FBA and community simulations

## 👷 Contributors

* Updated by Francisco Zorrilla in September 2024
* Expanded by Francisco Zorrilla and Arianna Basile in January 2024
* Originally developed by Arianna Basile and Kiran Patil in January 2023

## ☎️ Contact

Feel free to get in touch with us by raising an issue or sending an e-mail to Francisco Zorrilla (fz274@cam.ac.uk).
