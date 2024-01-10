# 2. Playing around with media conditions

## Authors
* Arianna Basile, MRC Toxicology Unit, University of Cambridge
* Francisco Zorrilla, MRC Toxicology Unit, University of Cambridge

## Learning outcomes

In this tutorial you will use [cobrapy](https://cobrapy.readthedocs.io/en/latest/) to learn how to:

- **2.1**: Modify growth medium of your reconstruction
- **2.2**: Perform gene essentiality analysis under different conditions
- **2.3**: Case study, simulate the Crabtree effect in yeast

## Setup


```python
# Import packages
import cobra

# Enable autocompleting with tab
%config Completer.use_jedi = False
```

## 2.1 Growth medium

The availability of nutrients has a major impact on metabolic fluxes and `cobrapy` provides some helpers to manage the exchanges between the external environment and your metabolic model. In experimental settings the “environment” is usually constituted by the growth medium, i.e. the concentrations of all metabolites and co-factors available to the modeled organism. However, constraint-based metabolic models only consider fluxes. Thus, you cannot simply use concentrations since fluxes have the unit `mmol / [gram of dry-cell weight * hour]` (i.e. concentration per gram dry weight of cells and hour).

Also, you are setting an upper bound for the particular import flux and not the flux itself. There are some crude approximations. For instance, if you supply 1 mmol of glucose every 24h to 1 gram of bacteria you might set the upper exchange flux for glucose to `1 mmol / [1 gDW * h]` since that is the nominal maximum that can be imported. There is no guarantee however that glucose will be consumed with that flux. Thus, the preferred data for exchange fluxes are direct flux measurements as the ones obtained from timecourse exa-metabolome measurements for instance.

So how does that look in COBRApy? The current growth medium of a model is managed by the medium attribute.


```python
# Import model
model_yeast=cobra.io.read_sbml_model("iMM904.xml.gz")

# Extract media information
medium=model_yeast.medium

# Show uptake bounds for extracellular media components
print("Default extracellular media components for yeast model: ")
medium
```

    Default extracellular media components for yeast model: 





    {'EX_fe2_e': 999999.0,
     'EX_glc__D_e': 10.0,
     'EX_h2o_e': 999999.0,
     'EX_h_e': 999999.0,
     'EX_k_e': 999999.0,
     'EX_na1_e': 999999.0,
     'EX_so4_e': 999999.0,
     'EX_nh4_e': 999999.0,
     'EX_o2_e': 2.0,
     'EX_pi_e': 999999.0}



This will return a dictionary that contains the upper flux bounds for all active exchange fluxes (the ones having non-zero flux bounds). Right now we see that we have enabled aerobic growth. 
Let's optimize to check the growth in the given medium.


```python
# Solve FBA
model_yeast.optimize()

# Show growth rate with glucose as carbon source
print("Growth rate in default medium: ",model_yeast.objective.value,"mmol/(gDW * h)")
```

    Growth rate in default medium:  0.28786570370401793 mmol/(gDW * h)


You can modify a growth medium of a model by assigning a dictionary to `model.medium` that maps exchange reactions to their respective upper import bounds. Let's edit the medium to have ethanol as carbon source instead of glucose.


```python
# Make a copy of the original model
model_yeastv2=model_yeast.copy()

# Extract media information
medium = model_yeastv2.medium

# Define new media uptake bounds
medium["EX_glc__D_e"] = 0
medium["EX_etoh_e"] = 1.0

# Update new media uptake bounds in model copy
model_yeastv2.medium = medium

# Show new uptake bounds
print("Modified extracellular media components for yeast model: ")
model_yeastv2.medium
```

    Modified extracellular media components for yeast model: 





    {'EX_etoh_e': 1.0,
     'EX_fe2_e': 999999.0,
     'EX_h2o_e': 999999.0,
     'EX_h_e': 999999.0,
     'EX_k_e': 999999.0,
     'EX_na1_e': 999999.0,
     'EX_so4_e': 999999.0,
     'EX_nh4_e': 999999.0,
     'EX_o2_e': 2.0,
     'EX_pi_e': 999999.0}



As we can see, the glucose in our medium has been replaced by ethanol. Let's check how it influences the growth rate.


```python
# Solve FBA with new uptake bounds
model_yeastv2.optimize()

# Show growth rate with ethanol as carbon source
print("Growth rate in modified medium with ethanol as carbon source: ",
      model_yeastv2.objective.value,"mmol/(gDW * h)")
```

    Growth rate in modified medium with ethanol as carbon source:  0.027542115993883887 mmol/(gDW * h)


The growth rate is reduced by around one order of magnitude when growing on ethanol compared to glucose. In conclusion, we have the same metabolic reconstruction (`iMM904`), simulated under two different media conditions using different carbon sources. `model_yeast` is growing with glucose as a carbon source while `model_yeastv2` is using ethanol.

### Questions

1. Can the yeast model grow using lactose as a carbon source? What about galactose? Was this an expected result?
2. What happens when you lower the uptake rate of ammonium? How about when it is set to 0? Was this expected?
3. Can you quantify the changes in growth rate induced by replacing ammonium with an amino acid of your choice to the media? e.g. `EX_val__L_e`

#### Tip #1
* Use the `model_yeast.exchanges` command to see which exchange reactions are present in the model.


```python
# Show first N exchange reactions
N=10
model_yeast.exchanges[:N]
```




    [<Reaction EX_epistest_SC_e at 0x7fead05d9d90>,
     <Reaction EX_epist_e at 0x7fead05d9c10>,
     <Reaction EX_ergst_e at 0x7fead093ea00>,
     <Reaction EX_ergstest_SC_e at 0x7fead0911e20>,
     <Reaction EX_13BDglcn_e at 0x7fead093e640>,
     <Reaction EX_etha_e at 0x7fead093eeb0>,
     <Reaction EX_2hb_e at 0x7fead093e9a0>,
     <Reaction EX_etoh_e at 0x7fead093eca0>,
     <Reaction EX_fe2_e at 0x7fead093ef70>,
     <Reaction EX_fecost_e at 0x7fead093edc0>]



#### Tip #2
* Use the [BiGG database](http://bigg.ucsd.edu/) to translate metabolite/reaction IDs. For example, click [here](http://bigg.ucsd.edu/search?query=EX_etha_e) to learn about the exchange reaction `EX_etha_e`.

### Answers

#### Question 1

Yeast cannot grow on lactose, as it cannot break it down into simple sugars [[ref](https://www.sciencedirect.com/science/article/pii/S2405471217303903)]. The model reflects this fact, as it is missing an exchange reaction for lactose. When we try to modify the uptake rate of lactose we get an error, because the model does not have an exchange reaction for lactose. On the other hand, yeast can indeed grow on galactose.


```python
# Make a copy of the original model
model_yeastv3=model_yeast.copy()

# Extract media information
medium = model_yeastv3.medium

# Define new media uptake bounds
medium["EX_glc__D_e"] = 0
medium["EX_etoh_e"] = 0
#medium["EX_lcts_e"] = 1
medium["EX_gal_e"] = 10

# Update new media uptake bounds in model copy
model_yeastv3.medium = medium

# Show new uptake bounds
print("Modified extracellular media components for yeast model with galactose as carbon source: ")
model_yeastv3.medium
```

    Modified extracellular media components for yeast model with galactose as carbon source: 





    {'EX_fe2_e': 999999.0,
     'EX_gal_e': 10,
     'EX_h2o_e': 999999.0,
     'EX_h_e': 999999.0,
     'EX_k_e': 999999.0,
     'EX_na1_e': 999999.0,
     'EX_so4_e': 999999.0,
     'EX_nh4_e': 999999.0,
     'EX_o2_e': 2.0,
     'EX_pi_e': 999999.0}




```python
# Solve FBA with new uptake bounds
model_yeastv3.optimize()

# Show growth rate with galactose as carbon source
print("\nGrowth rate in modified medium with galactose as carbon source:",
      model_yeastv3.objective.value,"mmol/(gDW * h)")
```

    
    Growth rate in modified medium with galactose as carbon source: 0.2878657037040199 mmol/(gDW * h)


#### Question 2

Lowering the ammonium uptake rate causes a decrease in growth rate. Setting the uptake rate to 0 results in near-zero growth rate. This is expected, as yeast needs a nitrogen source in order to grow [[ref](https://www.sciencedirect.com/science/article/pii/S2405471217303903)].


```python
# Make a copy of the original model
model_yeastv3=model_yeast.copy()

# Extract media information
medium = model_yeastv3.medium

# Define new media uptake bounds
medium["EX_glc__D_e"] = 10
medium["EX_etoh_e"] = 0
medium["EX_gal_e"] = 0
medium["EX_nh4_e"] = 1
#medium["EX_nh4_e"] = 0

# Update new media uptake bounds in model copy
model_yeastv3.medium = medium

# Show new uptake bounds
print("Modified extracellular media components for yeast model with less ammonium: ")
model_yeastv3.medium
```

    Modified extracellular media components for yeast model with less ammonium: 





    {'EX_fe2_e': 999999.0,
     'EX_glc__D_e': 10,
     'EX_h2o_e': 999999.0,
     'EX_h_e': 999999.0,
     'EX_k_e': 999999.0,
     'EX_na1_e': 999999.0,
     'EX_so4_e': 999999.0,
     'EX_nh4_e': 1,
     'EX_o2_e': 2.0,
     'EX_pi_e': 999999.0}




```python
# Solve FBA with new uptake bounds
model_yeastv3.optimize()

# Show growth rate with less nitrogen
print("Growth rate in modified medium with less ammonium:",
      model_yeastv3.objective.value,"mmol/(gDW * h)")
```

    Growth rate in modified medium with less ammonium: 0.17868438262760275 mmol/(gDW * h)


#### Question 3
Replacing ammonium with valine results in ~2.2% increase in growth rate.


```python
# Make a copy of the original model
model_yeastv3=model_yeast.copy()

# Extract media information
medium = model_yeastv3.medium

# Define new media uptake bounds
medium["EX_glc__D_e"] = 10
medium["EX_etoh_e"] = 0
medium["EX_gal_e"] = 0
medium["EX_nh4_e"] = 0
medium["EX_val__L_e"] = 1000

# Update new media uptake bounds in model copy
model_yeastv3.medium = medium

# Show new uptake bounds
print("Modified extracellular media components for yeast model with valine as nitrogen source: ")
model_yeastv3.medium
```

    Modified extracellular media components for yeast model with valine as nitrogen source: 





    {'EX_fe2_e': 999999.0,
     'EX_glc__D_e': 10,
     'EX_h2o_e': 999999.0,
     'EX_h_e': 999999.0,
     'EX_k_e': 999999.0,
     'EX_na1_e': 999999.0,
     'EX_so4_e': 999999.0,
     'EX_o2_e': 2.0,
     'EX_pi_e': 999999.0,
     'EX_val__L_e': 1000}




```python
# Solve FBA with new uptake bounds
model_yeastv3.optimize()

# Show growth rate with valine as nitrogen source
print("Growth rate in modified medium with valine as nitrogen source:",
      model_yeastv3.objective.value,"mmol/(gDW * h)")

# Calculate as percentage
print("Growth rate on valine as nitrogen source relative to ammonium:",
      (model_yeastv3.objective.value/model_yeast.objective.value)*100,"%")
```

    Growth rate in modified medium with valine as nitrogen source: 0.2942390011149862 mmol/(gDW * h)
    Growth rate on valine as nitrogen source relative to ammonium: 102.21398288471393 %


## 2.2 Gene essentiality

A gene is considered essential if restricting the flux of all reactions that depend on it to zero causes the objective, e.g., the growth rate, to also be zero, fall below a given threshold, or result in infeasible solutions.

Now we will use the function `find_essential_genes` to find essential genes in both conditions, extract the IDs of the genes, and compare the two lists. Can you do it on your own?

First lets check the model growing with ethanol as a carbon source:


```python
# Find essential genes for model growing under ethanol condition
essential_genes_v2=cobra.flux_analysis.variability.find_essential_genes(model_yeastv2)

# Get genes IDs
id_essential_genes_v2=[gene.id for gene in essential_genes_v2]

# Print number of esssential genes
print("Essential genes with ethanol carbon source:",len(id_essential_genes_v2))
```

    Essential genes with ethanol carbon source: 154


Next lets check the model growing with glucose as a carbon source:


```python
# Find essential genes for model growing under glucose condition
essential_genes=cobra.flux_analysis.variability.find_essential_genes(model_yeast)

# Get genes IDs
id_essential_genes=[gene.id for gene in essential_genes]

# Print number of essential genes
print("Essential genes with glucose carbon source:",len(id_essential_genes))
```

    Essential genes with glucose carbon source: 110


Check how large the intersect is between the two lists:


```python
a=list(set(id_essential_genes_v2) & set(id_essential_genes))
print("Intersect of essential genes between conditions:",len(a))
```

    Intersect of essential genes between conditions: 110


Lets check which genes are unique to each condition:


```python
print("Essential genes unique to ethanol growth condition:",
      len(list(set(id_essential_genes_v2) - set(id_essential_genes))))
list(set(id_essential_genes_v2) - set(id_essential_genes))
```

    Essential genes unique to ethanol growth condition: 44





    ['YKL060C',
     'YEL024W',
     'YDL181W',
     'YHR051W',
     'YBR039W',
     'YBR196C',
     'YLR295C',
     'YLR395C',
     'YLL041C',
     'YJR121W',
     'YLR038C',
     'YBL099W',
     'YGR183C',
     'YDR298C',
     'Q0130',
     'YPR191W',
     'YKL016C',
     'YER065C',
     'YKR097W',
     'YPL262W',
     'YDR529C',
     'YDL004W',
     'YDL067C',
     'YDR377W',
     'Q0250',
     'YGL187C',
     'YDR050C',
     'Q0275',
     'YBL045C',
     'YJL166W',
     'YCR012W',
     'Q0085',
     'Q0080',
     'YML081C_A',
     'YPL271W',
     'YLR377C',
     'YPL078C',
     'YMR256C',
     'YGL191W',
     'Q0045',
     'YHR001W_A',
     'YOR065W',
     'YFR033C',
     'Q0105']




```python
print("Essential genes unique to glucose growth condition:",
      len(list(set(id_essential_genes) - set(id_essential_genes_v2))))
list(set(id_essential_genes) - set(id_essential_genes_v2))
```

    Essential genes unique to glucose growth condition: 0





    []



### Questions

1. Why there are more essential genes for growth on ethanol than for growth on glucose? 
2. What pathways are these ethanol-growth essential genes associated with?
3. Choose an essential gene unique to growth in ethanol and search the literature for its function. Does it make sense to you that it is essential?

#### Tip 

For an overview of central carbon metabolism in yeast refer to the figure below, obtained from [this publication](https://journals.asm.org/doi/10.1128/mbio.02970-21).

![image.png](attachment:image.png)

### Answers

#### Question 1

Because growing on ethanol requires an extra set of enzymes which are not needed for growth on glucose. Growing on glucose requires enzymes for glycolysis only, while growing on ethanol additionally requires aerobic fermentation enzymes. 

#### Question 2

These genes are associated with the TCA cycle and oxidative phosphorylation (i.e. electron transport chain).

#### Question 3

A search for the gene [YLR395C](https://www.yeastgenome.org/locus/S000004387) reveals that it is indeed associated with the electron transport chain and cellular respiration.

## 2.3 Simulating the Crabtree effect 



Yeast central carbon metabolic pathways: all short-term Crabtree positive yeasts possess an upregulated aerobic (blue) and anaerobic (red) glycolytic pathway, even under fully aerobic conditions, when energy and carbon-source is limiting. At low glucose uptake rates, yeasts cells are purely respiring and there is a carbon-flux only through glycolysis (GF) and respiration (RF). Upon a sudden glucose excess condition (glucose pulse), the glycolytic flux will exceed the respiratory flux, which results in a fermentative flux (FF) and ethanol production. In other words, it appears as if slowly dividing short-term Crabtree positive cells, with little food around are already equipped with a strong energy producing apparatus, for rapid glucose consumption and energy production.

![image-3.png](attachment:image-3.png)
(Adapted from <a href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0116942">Hagman and Piškur, 2015</a>)

To simulate the Crabtree effect, we consider a fixed amount of glucose and different concentrations of oxygen in our medium. 


```python
# Create copy of model
model_yeastv3=model_yeast.copy()

# Define range of oxygen uptake rates
oxygen_uptakes=range(1,1000,10)
oxygen_uptakes = [ x / 100 for x in oxygen_uptakes]

# Define ethanol production vector
etoh_prod=[]

# Simulate under varying oxygen uptake rates, fixed glucose
for flux in oxygen_uptakes:
    medium = model_yeastv3.medium
    medium["EX_glc__D_e"] = 15
    medium["EX_o2_e"]=flux
    model_yeastv3.medium = medium
    opt=model_yeastv3.optimize()
    etoh_prod.append(opt.fluxes["EX_etoh_e"])

```

Plot results


```python
# Import plotting library
import matplotlib.pyplot as plt

# Define plot
x_axis=oxygen_uptakes
y_axis=etoh_prod
plt.plot(x_axis, y_axis)
plt.title('Crabtree effect in yeast model')
plt.xlabel('Oxygen uptake mmol / [gDW * h]')
plt.ylabel('Ethanol production')
plt.show()
```


    
![png](output_41_0.png)
    



```python
# Type your code here if you want to try similar analyses using different carbon sources, e.g. D-mannose

```

### Questions

1. Can you explain the observed behavior in the plot above?
2. Try playing around with other media components, for example, what happens if you keep oxygen fixed and vary glucose?
3. Literature search question: What are two similarities and two differences between the Crabtree and Warburg effects?

#### Tip

Check out [this publication](https://www.tandfonline.com/doi/full/10.1080/15384101.2018.1442622) describing the Crabtree and Warburg effects. In particular, have a look at Figure 1 shown below:

![image-2.png](attachment:image-2.png)

### Answers

#### Question 1

The plot shows the Crabtree effect. When glucose is present in excess and under anaerobic/low oxygen conditions, we see that ethanol production reaches a maximum. Furthermore, we see that as oxygen uptake increases, yeast reduces its fermentative flux but it does not reach zero.

#### Question 2

One of the conditions of the Crabtree effect is that the carbon source (i.e. glucose) has to be present in excess. In line with this precondition, we see that very little ethanol is produced (i.e. low fermentative flux) until glucose is no longer limiting.


```python
# Create copy of model
model_yeastv3=model_yeast.copy()

# Define range of oxygen uptake rates
glucose_uptakes=range(1,1000,10)
glucose_uptakes = [ x / 100 for x in glucose_uptakes]

# Define ethanol production vector
etoh_prod=[]

# Simulate under varying oxygen uptake rates, fixed glucose
for flux in glucose_uptakes:
    medium = model_yeastv3.medium
    medium["EX_o2_e"] = 10
    medium["EX_glc__D_e"]=flux
    model_yeastv3.medium = medium
    opt=model_yeastv3.optimize()
    etoh_prod.append(opt.fluxes["EX_etoh_e"])
    
# Define plot
x_axis=oxygen_uptakes
y_axis=etoh_prod
plt.plot(x_axis, y_axis)
plt.title('Crabtree effect in yeast model')
plt.xlabel('Glucose uptake mmol / [gDW * h]')
plt.ylabel('Ethanol production')
plt.show()
```

    /Users/fzorrilla/.local/lib/python3.8/site-packages/cobra/util/solver.py:554: UserWarning: Solver status is 'infeasible'.
      warn(f"Solver status is '{status}'.", UserWarning)



    
![png](output_44_1.png)
    


#### Question 3

Similarities:
* Both result in high fermentative flux even under aerobic conditions
* Both result in lower ATP yield

Differences:
* The Crabtree effect is observed in yeast cells, while the Warburg effect occurs in mammalian cells
* The Crabtree effect produces ethanol, while the Warburg effect produces lactate
