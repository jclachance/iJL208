###########################
### SOLVING FOR iJL208 ###
###########################
import json
import pandas as pd
import pulp
import itertools
import pdb
import re
import os
from tqdm import tqdm

def build_MIP_by_Cobrapy(model, growth_rate, essential_genes_file, parameters_file, regulator_genes_file, TU_Json_file, out_path='../data/minGenome', verbose=False, solver='CPLEX', iterations=10):
    M = 1000
    #Change variable names to comply former names
    me = model
    mu = growth_rate
    eg_f = essential_genes_file
    parameters_f = parameters_file
    reg_f = regulator_genes_file
    ############# sets ################################
    # TU
    with open(TU_Json_file) as data_file:    
        TUs = json.load(data_file)

    # essential genes        
    essential_genes = pd.read_csv(eg_f,index_col=0)
    essential_genes['gene'] = "u_G_" + essential_genes['gene'].astype(str)
    essential_genes = essential_genes['gene'].tolist()

    # regulator genes
    if reg_f != None:
        reg_genes = pd.read_csv(reg_f,index_col=0)
        reg_genes['gene'] = "u_G_" + reg_genes['gene'].astype(str)
        reg_genes = reg_genes['gene'].tolist()

    # ############# parameters ################################       
    df = pd.read_csv(parameters_f,index_col=0)

    test_all_genes = df["gene_or_promoter"].tolist()
    not_shared = []
    for gene in TUs.keys():
        if gene not in test_all_genes:
            not_shared.append(gene)

    df["gene_or_promoter"] = "u_G_" + df["gene_or_promoter"].astype(str)
    no_start = df[df['cannot_as_start']==1]["gene_or_promoter"].tolist()

    genes = df["gene_or_promoter"].tolist()

    end = df[['gene_or_promoter','start']].set_index('gene_or_promoter')\
                    .T.to_dict('list')
    start = df[['gene_or_promoter','start_if_select_as_start']]\
                    .set_index('gene_or_promoter').T.to_dict('list')

    reactions = [r_id.id for r_id in me.reactions]
    metabolites = [m_id.id for m_id in me.metabolites]

    ############# variables ################################
    v = pulp.LpVariable.dicts("v", reactions, 0, M, cat='Continuous')
    x = pulp.LpVariable.dicts("x", genes, cat='Binary')
    y = pulp.LpVariable.dicts("y", genes, cat='Binary')
    z = pulp.LpVariable.dicts("z", genes, cat='Binary') 
    # z can be defined as continuous

    ############# define model ################################
    lp_prob = pulp.LpProblem("MaxDeletion", pulp.LpMaximize)

    ############# objective ################################
    lp_prob += (pulp.lpSum([y[j]*end[j][0] for j in genes]) 
              - pulp.lpSum([x[j]*start[j][0] for j in genes])), "Max_length"

    def addReactionIndicator(lp_prob):
        for r in me.reactions:
            rgenes = r.genes
            GPR = r.gene_reaction_rule
            GPR = GPR.replace('\n','')
            GPR = GPR.replace('__10__','')

            if 's0001' in GPR: continue # not mapped gene in iJO1366
            if 'BG12900' in GPR: continue # not mapped gene in iYO844

            # pdb.set_trace()
            # no genes
            if len(rgenes) == 0:
                continue
            # single gene
            # if ('and' and 'AND' and 'or' and 'OR') not in GPR:
            if 'and' not in GPR \
            and 'AND' not in GPR \
            and 'or' not in GPR \
            and 'OR' not in GPR:
                # print GPR, genes
                assert(len(rgenes) == 1)
                for gene in rgenes:
                    gene = gene.id.replace('__10__','')
                    label = "knockout" + str(gene)
                    gene = "u_G_" + gene
                    lp_prob += v[r.id] - (1-z[gene])*M <= 0, \
                                        label + "_UB_" + r.id
                    lp_prob += v[r.id] - (1-z[gene])*(-M) >= 0, \
                                        label + "_LB_" + r.id

            # enzyme complex
            elif (('and' or 'AND') in GPR) and  (('or' or 'OR') not in GPR):
                # print genes
                # print GPR
                assert(len(rgenes) > 1)
                for gene in rgenes:
                    gene = gene.id.replace('__10__','').replace('(','').replace(')','')
                    label = "knockout_" + str(gene)
                    gene = "u_G_" + gene
                    lp_prob += v[r.id] - (1-z[gene])*M <= 0, \
                                        label + "_UB_" + r.id
                    lp_prob += v[r.id] - (1-z[gene])*(-M) >= 0, \
                                        label + "_LB" + r.id

            # isozymes
            elif (('and' or 'AND') not in GPR) and  (('or' or 'OR') in GPR):
                # print GPR
                lp_prob += v[r.id] - M <= 0, "knockout" + r.id + "Ori_UB"
                lp_prob += v[r.id] - (-M) >= 0, "knockout" + r.id + "Ori_LB"
                assert(len(rgenes) > 1)
                lp_prob += v[r.id] - M * pulp.lpSum(1-z['u_G_'+j.id.replace('__10__','')] \
                            for j in rgenes) <=0, "knockout" + r.id + '_UB'
                lp_prob += v[r.id] - (-M) * pulp.lpSum(1-z['u_G_'+j.id.replace('__10__','')] \
                            for j in rgenes) >=0, "knockout" + r.id + '_LB'
            # more complicated GPRs
            else:
#                 print r.id
#                 print GPR.split(' or ')                        
                proteins = [protein.replace("( ","").replace(" )","").split(' and ')\
                            for protein in GPR.split(' or ')]
    
                all_proteins = []        
                for protein in proteins:
                    mini = []
                    for prot in protein:
                        mini.append(prot.replace('(','').replace(')',''))
                    all_proteins.append(mini)
                    
                proteins = all_proteins
                commonGenes = set(proteins[0])
                # if len(gpr.proteins) > 1:
                for protein in proteins[1:]:
                    commonGenes.intersection_update(protein)
                nonCommonGenesList = []
                for protein in proteins:
                    nonCommonGenes = []
                    for gene in protein:
                        if gene not in commonGenes:
                            nonCommonGenes.append(gene)
                    nonCommonGenesList.append(nonCommonGenes)

                for gene in commonGenes:
                    # gene = gene.id
                    label = "knockout" + str(gene)
                    gene = "u_G_" + gene.replace('__10__','').replace('(','').replace(')','')
                    lp_prob += v[r.id] - (1-z[gene])*M <= 0, \
                                label + "_UB_" + r.id
                    lp_prob += v[r.id] - (1-z[gene])*(-M) >= 0, \
                                label + "_LB_" + r.id

                allCombination = list(itertools.product(*nonCommonGenesList))
                # print allCombination
                # print allCombination
                for i,genesC in enumerate(allCombination):
                    lp_prob += v[r.id] - M * pulp.lpSum(1-z['u_G_'+j.replace('__10__','')] \
                                for j in genesC) <=0,\
                        "knockout" + r.id + '_UB_' + str(i)
                    lp_prob += v[r.id] - (-M) * pulp.lpSum(1-z['u_G_'+j.replace('__10__','')] \
                                for j in genesC) >=0,\
                        "knockout" + r.id + '_LB_' + str(i)

    ############# constraints ################################
    if verbose:
        print("add reaction indicator")    
    addReactionIndicator(lp_prob)

    def get_S(model,mu):
        """build the stoichiometric matrix at a specific growth rate"""
            # intialize to 0
        # S = dok_matrix((len(self.metabolites), len(self.reactions)))
        S = {}
        # populate with stoichiometry
        for i, r in enumerate(model.reactions):
            for met, value in r._metabolites.items():
                #met_index = self.metabolites.index(met)
                if met.id not in S:
                    S[met.id] = {}
                if hasattr(value, "subs"):
                    S[met.id][r.id] = float(value.subs(mu, growth_rate))
                else:
                    S[met.id][r.id] = float(value)
        return S

    #### M-model constraints
    S = get_S(me, mu) # growth rate is 0.3
    # print S
    if verbose:
        print("add GSM constraint")
    # for i in metabolites:
    for i in S.keys():
        label = "mass_balance_%s"%i
        dot_S_v = pulp.lpSum([S[i][j] * v[j] for j in S[i].keys()])
        condition = dot_S_v == 0
        lp_prob += condition, label

    ###### cut in the genome
    if verbose:
        print("add cutting constraints")
    lp_prob += pulp.lpSum(y[j] for j in genes) == 1, "end"
    lp_prob += pulp.lpSum(x[j] for j in genes) == 1, "start"

    # cut genes between start and end
    # for i,gene in enumerate(genes):
    #     lp_prob += pulp.lpSum(x[j] for j in \
    #               genes[0:i+1]) - pulp.lpSum(y[j] \
    #                   for j in genes[0:i+1]) - z[gene] == 0,\
    #                        'indicator' + str(gene)

    # A = pulp.LpAffineExpression()

    # for i,gene in enumerate(genes):
    #     A.addterm(x[gene],1) #pulp.lpSum(x[j] for j in genes[0:i+1])
    #     A.addterm(y[gene],-1)
    #     lp_prob += A - z[gene] == 0,'indicator' + str(gene)

    lp_prob += x[genes[0]] - y[genes[0]] == z[genes[0]], 'indicator' + str(genes[0])
    for i,gene in enumerate(genes):
        if i == 0: continue
        lp_prob += z[genes[i-1]] + x[gene] - y[gene] == z[gene],'indicator' + str(gene)

    ##### TUs 
    if verbose:
        print("add TU constraint")
    for gene,promoters in TUs.items():
        if gene in not_shared: continue
        gene = 'u_G_' + gene
        len_pro = len(promoters)
        #print(gene, promoters)
        lp_prob += z[gene] - pulp.lpSum(z['u_G_'+j] for j in promoters) + \
                    (len_pro - 1) >= 0,'TU_all_'+gene
        for pro in promoters:
            pro = 'u_G_' + pro
            lp_prob += z[gene] - z[pro] <=0, 'TU_'+gene+'_'+pro

    ##### some overlapped region cannot be selected as the start of deletion
    if verbose:
        print("add no start and essential genes")
    for gene in no_start:
        lp_prob += x[gene] == 0, 'no_start_'+gene

    # knock out transcription of cutted genes
    for gene in genes:
        label = "knockout" + str(gene)
        # pdb.set_trace()
        if gene in v.keys():
            lp_prob += v[gene] - (1-z[gene])*M <= 0, label

    ##### add essential genes that cannot be deleted
    for eg in essential_genes:
        if eg in genes:
            lp_prob += z[eg] == 0

    ##### add regulation genes that cannot be deleted
    if reg_f != None:
        for eg in reg_genes:
            # remove the part joint with essential genes
            if (eg in genes) and (eg not in essential_genes):
                lp_prob += z[eg] == 0

    ##### reaction bounds
    for r in me.reactions:
        # (lb,up) = me.bounds[r_id]
        v[r.id].lowBound = r.lower_bound
        v[r.id].upBound = r.upper_bound
    
    v['BIOMASS_step3_c'].lowBound = mu

    v['BIOMASS'].lowBound = 0
    v['BIOMASS'].upBound = 0
    
    v['BIOMASS_step1_c'].lowBound = 0
    v['BIOMASS_step1_c'].upBound = 0
    
    v['BIOMASS_step2_c'].lowBound = 0
    v['BIOMASS_step2_c'].upBound = 0

    # lp file is somtime too larget to write
    # lp_prob.writeLP(lpfilename)

    # orignial implementation in the paper was calling cplex from C++ directly
    # call eternal compled cpp excutable to solve is a better option
    # it is implemented in codebase/mingenome_ecoli.cpp

    # current test version of using python to call the optimization
    # options = [epgap, epagap, epint, epopt, eprhs]
    if solver == 'gurobi':
        GUROBI_CMD_OPTIONS = [('Threads', 8), ('TimeLimit', 1800), ('FeasibilityTol',1E-9),
                        ('OptimalityTol',1E-9),('IntFeasTol',1E-9),
                        ('MIPGapAbs', 0), ('MIPGap', 0), ('CliqueCuts', 2)]
        pulp_solver = pulp.solvers.GUROBI_CMD(path=None, keepFiles=0, mip=1, msg=0,
                                options=GUROBI_CMD_OPTIONS)
    elif solver == 'CPLEX':
        pulp_solver = pulp.solvers.CPLEX(path=None, keepFiles=0, mip=1,\
            msg=1, options=['mip tolerances mipgap 0', \
            'mip tolerances absmipgap 0', 'mip tolerances integrality 0',\
            'simplex tolerances optimality 1E-9',\
            'simplex tolerances feasibility 1E-9',], timelimit=1200)
    elif solver == 'GLPK':
        pulp_solver = pulp.solvers.GLPK(path=None, keepFiles=0, mip=1,\
            msg=1, options=['mip tolerances mipgap 0', \
            'mip tolerances absmipgap 0', 'mip tolerances integrality 0',\
            'simplex tolerances optimality 1E-9',\
            'simplex tolerances feasibility 1E-9',])
    else:
        raise ValueError('Solver name not compatible')

    x_list = []
    y_list = []
    status = []

    def iterate_solve(lp_prob,iter_count):
        lp_prob.solve(pulp_solver)
        if verbose:
            print("----------- " + str(iter_count) + " ------------")
        status.append(pulp.LpStatus[lp_prob.status])
        if verbose:
            print("Status:", pulp.LpStatus[lp_prob.status])
        for v in lp_prob.variables():
            if "x_u_G_" in v.name and v.varValue == 1:               
                xname = v.name.replace("x_","")
                xname = xname.replace('_','-')
                xname = xname.replace("PM-","PM_")
                xname = xname.replace('u-','u_')
                xname = xname.replace('G-','G_')
                #print(xname,v.name)
                lp_prob += x[xname] == 1
                if xname not in x_list: 
                    x_list.append(xname)
            if "y_u_G_" in v.name and v.varValue == 1:
                yname = v.name.replace("y_","")
                yname = yname.replace('_','-')
                yname = yname.replace("PM-","PM_")
                yname = yname.replace('u-','u_')
                yname = yname.replace('G-','G_')
                lp_prob += y[yname] == 1
                if yname not in y_list: 
                    y_list.append(yname)
        rhs = iter_count + 1
        lp_prob.constraints['start'].changeRHS(rhs)
        lp_prob.constraints['end'].changeRHS(rhs)
        return lp_prob

    for iter_count in range(1,iterations):
        #Updates the lp_prob at each iteration
        lp_prob = iterate_solve(lp_prob,iter_count)
    
    #Write the final results  
    out_file = 'deletion_results_' + str(iteration-1) + '.csv'
    writing_path = os.path.join(out_path, out_file)
    pd.DataFrame({'start': x_list, 'end':y_list, 'status':status}).to_csv(writing_path)

#### analyze result
def get_all_deletions(result_df, genes_and_promoters):
    #Get the start and end location, and the span of them
    all_deletions = []
    for i, row in result_df.iterrows():
        start_element = row['start'].replace('u_G_','')
        end_element = row['end'].replace('u_G_','')
        #Find start and end in genome
        for j, line in genes_and_promoters.iterrows():
            if start_element == line['gene_or_promoter']:
                start = line['start']
            if end_element == line['gene_or_promoter']:
                end = line['end']
        all_deletions.append((start,end, abs(start-end)))
    deletions_loc = pd.DataFrame.from_records(all_deletions, columns=['start_loc','end_loc','length'])
    return all_deletions
    
def get_genes_in_results(all_deletions, genes_and_promoters):
    #Get all the genes in the results
    deleted_genes = []
    for t in all_deletions:
        # '+' strand deletion
        if t[1] - t[0] > 0:
            start = t[0]
            end = t[1]
        # '-' strand deletions
        elif t[1] - t[0] < 0:
            start = t[1]
            end = t[0]
        #Find the genes within those boundaries
        deleted_genes.append([g for g in genes_and_promoters['gene_or_promoter'][(genes_and_promoters['start'] > start)\
                                                            & (genes_and_promoters['end'] < end)\
                                                            & (genes_and_promoters['class'] == 'genes')]])
    all_deleted_genes = []
    for l in deleted_genes:
        for g in l:
            all_deleted_genes.append(g)
    return list(set(all_deleted_genes))
    
def calculate_mcc(all_deleted_genes, comparison_syn3):
    from math import sqrt
    def get_confusion_matrix(all_deleted_genes, new_baby_sheet):
        #Make the comparisons now 
        #Number of deleted genes absent from syn3.0 (true positives)
        true_positives = set(all_deleted_genes).intersection(set(new_baby_sheet['locus_tag'][new_baby_sheet['syn3.0'] == 'thrash'].to_list()))
        #Number of deleted genes that are in syn3.0 (false positives)
        false_positives = set(all_deleted_genes).intersection(set(new_baby_sheet['locus_tag'][new_baby_sheet['syn3.0'] == 'keep'].to_list()))
        #Number of non-deleted genes that are in syn3.0 (true negatives)
        all_florum_genes = set(new_baby_sheet['locus_tag'].to_list())
        non_deleted_genes = all_florum_genes.difference(set(all_deleted_genes))
        true_negatives = non_deleted_genes.intersection(set(new_baby_sheet['locus_tag'][new_baby_sheet['syn3.0']=='keep'].to_list()))
        #Number of non-deleted genes that are missing in syn3.0 (false negatives)
        false_negatives = non_deleted_genes.intersection(set(new_baby_sheet['locus_tag'][new_baby_sheet['syn3.0']=='thrash']))
        
        return len(true_positives), len(false_positives), len(true_negatives), len(false_negatives)
    
    tp, fp, tn, fn = get_confusion_matrix(all_deleted_genes, comparison_syn3)
    num = float((tp*tn)-(fp*fn))
    denom = float(sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
    mcc = num/denom
    return mcc

def get_deletion_results(max_deletion_df, genes_and_promoters, comparison_syn3):
    all_deletion_results, old_all_deleted_genes = [], []
    all_deletions = get_all_deletions(max_deletion_df, genes_and_promoters)
    for i in tqdm(range(len(max_deletion_df))):
        result_df = max_deletion_df.iloc[:i,:]
        if result_df.empty:
            pass
        else:
            new_all_deleted_genes = get_genes_in_results(all_deletions[:i], genes_and_promoters)
            deleted_genes_in_deletion = list(set(new_all_deleted_genes).difference(set(old_all_deleted_genes)))
            old_all_deleted_genes = new_all_deleted_genes
            mcc = calculate_mcc(old_all_deleted_genes, comparison_syn3)
            # Generate the deletions results for this iteration
            all_deletion_results.append((len(old_all_deleted_genes), 
                                         deleted_genes_in_deletion,
                                         sum([t[2] for t in all_deletions[:i]]),
                                         mcc))
    
    return all_deletion_results
