import cobra
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
from cobra.flux_analysis import pfba
from cobra.flux_analysis import single_gene_deletion

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sb

import numpy as np
import pandas as pd
import os
import glob

#Utilitary functions
def save_model(model):
    import time
    timestr = time.strftime("%Y_%m_%d_%H:%M:%S")
    path = '/ip29/jclachan/Doctorat/iJL207/Model_versions/'
    model_name = 'iJL' + str(len(model.genes)) + '_' + timestr +'.json'
    cobra.io.save_json_model(model,os.path.join(path,model_name))
    
    
def load_latest_model():
    import glob
    #Path to where all the models are stored
    path = '../Model_versions/*'
    #All files in that path
    list_of_files = glob.glob(path) # * means all if need specific format then *.csv
    #Get the latest
    latest_file = max(list_of_files, key=os.path.getctime)
    return cobra.io.load_json_model(latest_file)

def make_metabolite(metab_id,metab_formula,metab_name,compartment):
    from cobra import Metabolite
    metab = Metabolite(metab_id,
                       formula=metab_formula,
                       name=metab_name,
                       compartment=compartment)
    
    return metab

def make_reaction(rxn_id,rxn_name,rxn_subsystem,
                  rxn_lower_bound = -1000., rxn_upper_bound =1000):
    from cobra import Reaction
    rxn = Reaction(rxn_id)
    rxn.name = rxn_name
    rxn.susbsystem = rxn_subsystem
    rxn.lower_bound = rxn_lower_bound
    rxn.upper_bound = rxn_upper_bound
    
    return rxn

from cobra.util import linear_reaction_coefficients
from BOFdat.util.update import _get_biomass_objective_function

def _assess_metab_solvability(m, model):
    # Identify the list of metabolites that do not prevent the model to solve when added to the BOF
    model.solver = 'gurobi'
    biomass = _get_biomass_objective_function(model)
    biomass.remove_from_model()
    BIOMASS = Reaction('BIOMASS')
    model.add_reactions([BIOMASS])
    model.reactions.BIOMASS.add_metabolites({m: -1.})
    model.reactions.BIOMASS.objective_coefficient = 1.
    solution = model.slim_optimize()
    # If the model can produce that metabolite
    if solution > 1e-9:
        return (m,True)
    else:
        return (m,False)
    
def _assess_rxn_solvability(r, model):
    # Identify the list of metabolites that do not prevent the model to solve when added to the BOF
    model.solver = 'gurobi'
    try:
        biomass = _get_biomass_objective_function(model)
        biomass.objective_coefficient = 0.
    except: 
        pass
    model.reactions.get_by_id(r.id).objective_coefficient = 1.
    solution = model.slim_optimize()
    # If the model can produce that metabolite
    if solution > 1e-9:
        return (r.id,True)
    else:
        return (r.id,False)
    
def set_media(model, media):
    """
    Media is a list of reaction identifiers
    """
    CONSTRAINTS = ['EX_ac_e', 'EX_lac__L_e', 'EX_o2_e']
    #Set all bounds to 0 without changing the constraints
    for r in model.exchanges:
        if r.id not in CONSTRAINTS:
            r.lower_bound = 0.
        if r.id in media:
            r.lower_bound = -10.5

    return model

def change_energy_source(rxn_id, model, energy_exchanges):
    for r in model.reactions:
        if r.id == rxn_id:
            print(rxn_id)
            r.lower_bound =-18.;r.upper_bound=1000.
        elif r.id in energy_exchanges:
            r.lower_bound = 0.;r.upper_bound = 1000.
        else:
            r.lower_bound = -10.;r.upper_bound = 1000.
        #Set lactate and acetate secretion as constraints
        if r.id == 'EX_lac__L_e':
            r.lower_bound = 7.0; r.upper_bound = 12.
        if r.id == 'EX_ac_e':
            r.lower_bound = 2.0; r.upper_bound = 1000.

    return model
