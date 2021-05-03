#! /usr/bin/python3
#-*-coding: utf-8-*-

import cobra
from cobra import Model, Reaction, Metabolite
import tempfile
from pprint import pprint
from cobra.io import write_sbml_model, validate_sbml_model

import cobra.test
import os
from os.path import join

from cobra.util.solver import linear_reaction_coefficients

from cobra.flux_analysis import flux_variability_analysis



import pandas
from time import time

import cobra.test
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

from cobra.flux_analysis import production_envelope

from cobra.test import create_test_model
from cobra.sampling import sample
from cobra.sampling import OptGPSampler, ACHRSampler

import numpy as np 



cobra_config = cobra.Configuration()

cobra_config.lower_bound
cobra_config.upper_bound
cobra_config.bounds
cobra_config.bounds = -10, 20
cobra.Reaction("R1")
cobra.Reaction("R2", lower_bound=None)



#######################################################################################################################################################
################################################################# Building a model ####################################################################
#######################################################################################################################################################



############################################################################
####################### Model, reaction and metabolites ####################
############################################################################



############################################################
### 1st step : Create the model and reaction
############################################################

model = Model('example_model')
reaction = Reaction('R_30AS140')
reaction.name = '3 oxoacyl acyl carrier bound protein synthase n C140'
reaction.subsystem = 'Cell Envelope Biosynthesis'
reaction.lower_bound = 0 # default
reaction.upper_bound = 1000 # default

############################################################
### 2nd step : Create the metabolites
############################################################

ACP_c = Metabolite(
    'ACP_c',
    formula = 'C11H21N2O7PRS',
    name = 'acyl-carrier-protein',
    compartment = 'c')
omrsACP_c = Metabolite(
    'M3omrsACP_c',
    formula = 'C25H45N2O9PRS',
    name = '3-Oxotetradecanoyl-acyl-carrier-protein',
    compartment = 'c'
)
co2_c = Metabolite('co2_c', formula = 'CO2', name = 'CO2', compartment = 'c')
malACP_c = Metabolite(
    'malACP_c',
    formula = 'C14H22N2O10PRS',
    name = 'Malonyl-acyl-carrier-protein',
    compartment = 'c'
)
h_c = Metabolite('h_c', formula = 'H', name = 'H', compartment = 'c')
ddcaACP_c = Metabolite(
    'ddcaACP_c',
    formula = 'C23H43N2O8PRS',
    name = 'Dodecanoyl-ACP-n-C120ACP',
    compartment = 'c'
)

# Side note (SId) :
# highly recommended for metabolites, reactions and genes => "make sure" that they are valid SBML identifiers (SId)
# SId = data type derived from the basic XML typestring, but with restrictions about the characters permitted and the seq in which those characters may 
# appear
# Main restriction : these id cannot start with a number
# SId => serialization to SBML


########################################################
### 3rd step : Adding metabolites to a reaction
########################################################

# => Uses a dictionary : metabolite as key and stoichiometric coefficient as value

reaction.add_metabolites({
    malACP_c : -1.0,
    h_c: -1.0,
    ddcaACP_c : -1.0,
    co2_c : 1.0,
    ACP_c : 1.0,
    omrsACP_c: 1.0
})
print(reaction.reaction) # gives a string representation of the reaction

reaction.gene_reaction_rule = '( STM2378 or STM1197 )' 
# gene_reaction_rule = boolean representation of the gene requirements for the reaction to be active
print(reaction.genes)

### Model still empty for now
print(f"{len(model.reactions)} reactions initially")
print(f"{len(model.metabolites)} metabolites initially")
print(f"{len(model.genes)} genes initially")


########################################################
### 4th step : Add the reaction to the model
########################################################

model.add_reactions([reaction]) # objects added to the model

print(f'{len(model.reactions)} reactions')
print(f'{len(model.metabolites)} metabolites')
print(f'{len(model.genes)} genes')

# Iterate through the model objects to observe the contents
print("Reactions")
print("---------")
for x in model.reactions:
    print("%s : %s" % (x.id, x.reaction)) 

# The %s token allows us to insert (and potentially format) a string. Notice that the %s token is replaced by whatever we pass to the string after 
# the % symbol. Notice also that we are using a tuple here as well (when you only have one string using a tuple is optional) to illustrate that multiple 
# strings can be inserted and formatted in one statement.

print("")
print("Metabolites")
print("-----------")
for x in model.metabolites:
    print('%9s : %s' % (x.id, x.formula))

print("")
print("Genes")
print("-----")
for x in model.genes:
    associated_ids = (i.id for i in x.reactions)
    print("%s is associated with reactions: %s" %
          (x.id, "{" + ", ".join(associated_ids) + "}"))

# 9 is field width. If you pass string as Hi, this will print Hi and seven spaces.
# If you pass string as Hai, this will print Hai and six spaces.



###########################################################################
################################ Objective ################################
###########################################################################



# Set the objective of the model
model.objective = 'R_30AS140'
print(model.objective.expression)
print(model.objective.direction)



###########################################################################
############################# Model validation ############################
###########################################################################



# Validate and export the model to SBML for exchange with other tools
with tempfile.NamedTemporaryFile(suffix=' .xml') as f_sbml:
    write_sbml_model(model, filename = f_sbml.name)
    report = validate_sbml_model(filename = f_sbml.name)

pprint(report)
# lists are empty => no COBRA or SBML errors or warnings => the model seems to be valid



###########################################################################
######################### Exchanges, sinks and demands ####################
###########################################################################



# pre-defined boundary reactions; defined on metabolites
print("exchanges", model.exchanges)
print("demands", model.demands)
print("sinks", model.sinks)

model.add_metabolites([
    Metabolite(
        'glycogen_c',
        name='glycogen',
        compartment='c'
    ),
    Metabolite(
        'co2_e',
        name='CO2',
        compartment='e'
    ),
])

# create exchange reaction
print(model.add_boundary(model.metabolites.get_by_id("co2_e"), type="exchange"))

# create exchange reaction
print(model.add_boundary(model.metabolites.get_by_id("glycogen_c"), type="sink")) #write type="demand" instead of type="sink" to create a demand reaction

# Now we have an additional exchange and sink reaction in the model
print("exchanges", model.exchanges)
print("sinks", model.sinks)
print("demands", model.demands)

# boundary reactions
print(model.boundary) # to get the information about all boundary reactions => property boundary

# metabolic reactions
print(set(model.reactions) - set(model.boundary))



#########################################################################################################################################################
########################################################### Reading and Writing Models ##################################################################
#########################################################################################################################################################



# in SBML
data_dir = cobra.test.data_dir

print("mini test files: ")
print(", ".join(i for i in os.listdir(data_dir) if i.startswith("mini")))

textbook_model = cobra.test.create_test_model("textbook")
#ecoli_model = cobra.test.create_test_model("ecoli")
#salmonella_model = cobra.test.create_test_model("salmonella")



#############################################################################
################################## SBML #####################################
#############################################################################



# SBML = System Biology Markup Langage
# = XML-based standard format for distributing models which has support for COBRA models through the FBC extension version 2 (flux balance constraints)
# WARNING !!! SId requirements
print(cobra.io.read_sbml_model(join(data_dir, "mini_fbc2.xml"))) # does not generate warning
print(cobra.io.write_sbml_model(textbook_model, "test_fbc2.xml")) # does not generate warning, and generates a file, named test_fbc2.xml

#print(cobra.io.read_sbml_model(join(data_dir, "mini_cobra.xml"))) 

print(cobra.io.write_sbml_model(textbook_model, "test_cobra.xml"))



#############################################################################
################################### JSON ####################################
#############################################################################



print(cobra.io.load_json_model(join(data_dir, "mini.json")))
print(cobra.io.save_json_model(textbook_model, "test.json"))



#############################################################################
################################### YAML ####################################
#############################################################################



print(cobra.io.load_yaml_model(join(data_dir, "mini.yml")))
print(cobra.io.save_yaml_model(textbook_model, "test.yml"))



#############################################################################
################################## MATLAB ###################################
#############################################################################



print(cobra.io.load_matlab_model(join(data_dir, "mini.mat"), variable_name="mini_textbook"))
print(cobra.io.load_matlab_model(join(data_dir, "mini.mat")))
print(cobra.io.save_matlab_model(textbook_model, "test.mat"))



#########################################################################################################################################################
################################################################# Simulating with FBA ###################################################################
#########################################################################################################################################################



model = cobra.test.create_test_model("textbook")



###########################################################################
############################# Running with FBA ############################
###########################################################################



solution = model.optimize() # return a solution object
print(solution)
print(solution.objective_value)

#%%time # %%time prints the wall time for the entire cell whereas %time gives you the time for first line only. => part of IPython
#print(model.optimize().objective_value)

#%%time
#print(model.slim_optimize())


#######################################
### 1st step : Analyzing FBA solution
#######################################

print(model.summary())
print(model.metabolites.nadh_c.summary())
print(model.metabolites.atp_c.summary())


#######################################
### 2nd step : Changing the objectives
#######################################


biomass_rxn = model.reactions.get_by_id("Biomass_Ecoli_core")
print(linear_reaction_coefficients(model))

# change the objective to ATPM
model.objective = "ATPM"

# The upper bound should be 1000, so that we get the actual optimal value
model.reactions.get_by_id("ATPM").upper_bound = 1000.
print(linear_reaction_coefficients(model))
print(model.optimize().objective_value)


#######################################
### 3rd step : Running FVA
#######################################


# Flux variability analysis => finds the ranges of each metabolite flux at the optimum
#print(flux_variability_analysis(model, model.reactions[:10]))
# Setting parameter fraction_of_optimium=0.90 would give the flux ranges for reactions at 90% optimality.
#print(cobra.flux_analysis.flux_variability_analysis(model, model.reactions[:10], fraction_of_optimum=0.9))
loop_reactions = [model.reactions.FRD7, model.reactions.SUCDi]
#print(flux_variability_analysis(model, reaction_list=loop_reactions, loopless=False))
#print(flux_variability_analysis(model, reaction_list=loop_reactions, loopless=True))


### Running FVA in summary methods
model.optimize()
#print(model.summary(fva=0.95))
#print(model.metabolites.pyr_c.summary(fva=0.95))



########################################################################
########################### Running pFBA ###############################
########################################################################



# Parsimonious FBA => finds a flux distribution which gives the optimal growth rate, but minimizes the total  sum of flux
# => solve 2 sequential linear programs
model.objective = 'Biomass_Ecoli_core'
fba_solution = model.optimize()
pfba_solution = cobra.flux_analysis.pfba(model)
print(abs(fba_solution.fluxes["Biomass_Ecoli_core"] - pfba_solution.fluxes["Biomass_Ecoli_core"]))



########################################################################
########################## Running geometric FBA #######################
########################################################################



# finds a unique optimal flux distribution 
geometric_fba_sol = cobra.flux_analysis.geometric_fba(model)
print(geometric_fba_sol)



#########################################################################################################################################################
################################################################# Simulating deletions ##################################################################
#########################################################################################################################################################



cobra_model = cobra.test.create_test_model("textbook")
ecoli_model = cobra.test.create_test_model("ecoli")



###############################################################################
### Knocking-out single genes and reactions
###############################################################################



print('complete model: ', cobra_model.optimize())
with cobra_model:
    cobra_model.reactions.PFK.knock_out()
    print('pfk knocked out: ', cobra_model.optimize())

print('complete model: ', cobra_model.optimize())
with cobra_model:
    cobra_model.genes.b1723.knock_out()
    print('pfkA knocked out: ', cobra_model.optimize())
    cobra_model.genes.b3916.knock_out()
    print('pfkB knocked out: ', cobra_model.optimize())



################################################################################
### Single deletions
################################################################################



deletion_results = single_gene_deletion(cobra_model)
#print(single_gene_deletion(cobra_model, cobra_model.genes[:20]))
#print(single_reaction_deletion(cobra_model, cobra_model.reactions[:20]))



#################################################################################
### Double deletions
#################################################################################



#print(double_gene_deletion(cobra_model, cobra_model.genes[-5:]).round(4))

start = time()  # start timer()
double_gene_deletion(
    ecoli_model, ecoli_model.genes[:25], processes=2)
t1 = time() - start
print("Double gene deletions for 200 genes completed in "
      "%.2f sec with 2 cores" % t1)

start = time()  # start timer()
double_gene_deletion(
    ecoli_model, ecoli_model.genes[:25], processes=1)
t2 = time() - start
print("Double gene deletions for 200 genes completed in "
      "%.2f sec with 1 core" % t2)

print("Speedup of %.2fx" % (t2 / t1)) # f for floating point numbers; %.2f to display only 2 decimals; x for octal hexadecimal numbers


print(double_reaction_deletion(cobra_model, cobra_model.reactions[2:7]).round(4))



###############################################################################
### Accessing individual deletion results
###############################################################################



single = single_reaction_deletion(cobra_model)
double = double_reaction_deletion(cobra_model)

print(single.knockout["ATPM"])
print(double.knockout[{"ATPM", "TKT1"}])

atpm = cobra_model.reactions.ATPM
tkt1 = cobra_model.reactions.TKT1
pfk = cobra_model.reactions.PFK

print(single.knockout[atpm, tkt1, pfk])
print(double.knockout[{atpm, tkt1}, {atpm, pfk}, {atpm}])



###############################################################################
### Production envelopes
###############################################################################



# => show distinct phases of optimal growth with different use of two different substrates
model = cobra.test.create_test_model("textbook")
prod_env = production_envelope(model, ["EX_glc__D_e", "EX_o2_e"])
#print(prod_env.head())
prod_env = production_envelope(model, ["EX_o2_e"], objective="EX_ac_e", carbon_sources="EX_glc__D_e")
print(prod_env.head())

# In jupyter

%matplotlib inline
prod_env.plot(kind='line', x='EX_o2_e', y='carbon_yield_maximum');



#########################################################################################################################################################
################################################################# Flux sampling #########################################################################
#########################################################################################################################################################



################################################################################
### Basic usage
################################################################################



model = create_test_model("textbook")
s = sample(model, 100)
#print(s.head())

print("One process:")
%time s = sample(model, 1000) # synthax error (%)
print("Two processes:")
%time s = sample(model, 1000, processes=2)


s = sample(model, 100, method="achr")



################################################################################
### Advanced usage
################################################################################



achr = ACHRSampler(model, thinning=10)
optgp = OptGPSampler(model, processes=4)
#print(achr)
#print(optgp)


################################################################################
### Sampling and validation
################################################################################



# => function sample()
s1 = achr.sample(100)
s2 = optgp.sample(100)
#print(s1)
#print(s2)

bad = np.random.uniform(-1000, 1000, size=len(model.reactions))
#print(achr.validate(np.atleast_2d(bad)))
#print(achr.validate(s1))
s1_valid = s1[achr.validate(s1) == "v"]
#print(len(s1_valid))



################################################################################
### Batch sampling
################################################################################



counts = [np.mean(s.Biomass_Ecoli_core > 0.1) for s in optgp.batch(100, 10)]
#print("Usually {:.2f}% +- {:.2f}% grow...".format(np.mean(counts) * 100.0, np.std(counts) * 100.0))



################################################################################
### Adding constraints
################################################################################



co = model.problem.Constraint(model.reactions.Biomass_Ecoli_core.flux_expression, lb=0.1)
model.add_cons_vars([co])
s = sample(model, 10)
print(s.Biomass_Ecoli_core)



#########################################################################################################################################################
################################################################### Loopless FBA ########################################################################
#########################################################################################################################################################



#%matplotlib inline
#import plot_helper

import cobra.test
from cobra import Reaction, Metabolite, Model
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
from cobra.flux_analysis import pfba
import pandas



###########################################################################
### Loopless solution
###########################################################################


salmonella = cobra.test.create_test_model('salmonella')
nominal = salmonella.optimize()
loopless = loopless_solution(salmonella)
df = pandas.DataFrame(dict(loopless=loopless.fluxes, nominal=nominal.fluxes))
#df.plot.scatter(x='loopless', y='nominal')



###########################################################################
### Loopless model
###########################################################################



#plot_helper.plot_loop() => import matplotlib.pyplot

model = Model()
model.add_metabolites([Metabolite(i) for i in "ABC"])
model.add_reactions([Reaction(i) for i in ["EX_A", "DM_C", "v1", "v2", "v3"]])

model.reactions.EX_A.add_metabolites({"A": 1})
model.reactions.DM_C.add_metabolites({"C": -1})

model.reactions.v1.add_metabolites({"A": -1, "B": 1})
model.reactions.v2.add_metabolites({"B": -1, "C": 1})
model.reactions.v3.add_metabolites({"C": -1, "A": 1})

model.objective = 'DM_C'


with model:
    add_loopless(model)
    solution = model.optimize()
print("loopless solution: status = " + solution.status)
print("loopless solution flux: v3 = %.1f" % solution.fluxes["v3"])

solution = pfba(model)
print("parsimonious solution: status = " + solution.status)
print("loopless solution flux: v3 = %.1f" % solution.fluxes["v3"])

model.reactions.v3.lower_bound = 1
with model:
    add_loopless(model)
    try:
        solution = model.optimize()
    except:
        print('model is infeasible')

solution = pfba(model)
print("parsimonious solution: status = " + solution.status)
print("loopless solution flux: v3 = %.1f" % solution.fluxes["v3"])



############################################################################
### Method
############################################################################



#########################################################################################################################################################
############################################################## Consistency testing ######################################################################
#########################################################################################################################################################



import cobra
test_model = cobra.Model("test_model")
v1 = cobra.Reaction("v1")
v2 = cobra.Reaction("v2")
v3 = cobra.Reaction("v3")
v4 = cobra.Reaction("v4")
v5 = cobra.Reaction("v5")
v6 = cobra.Reaction("v6")

test_model.add_reactions([v1, v2, v3, v4, v5, v6])

v1.reaction = "-> 2 A"
v2.reaction = "A <-> B"
v3.reaction = "A -> D"
v4.reaction = "A -> C"
v5.reaction = "C -> D"
v6.reaction = "D ->"

v1.bounds = (0.0, 3.0)
v2.bounds = (-3.0, 3.0)
v3.bounds = (0.0, 3.0)
v4.bounds = (0.0, 3.0)
v5.bounds = (0.0, 3.0)
v6.bounds = (0.0, 3.0)

test_model.objective = v6


##########################################################################
### Using FVA
##########################################################################


#print(cobra.flux_analysis.find_blocked_reactions(test_model))


##########################################################################
### Using fastcc
##########################################################################


consistent_model = cobra.flux_analysis.fastcc(test_model)
print(consistent_model.reactions)
