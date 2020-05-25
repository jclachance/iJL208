###########################
### SOLVING FOR iJL208 ###
###########################

def build_MIP_by_Cobrapy(me,mu,
    eg_f,
    parameters_f,
    reg_f,
    TU_Json_file):

    M = 1000

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
    print("add reaction indicator")    
    addReactionIndicator(lp_prob)

    #### ME model constraint
    S = get_S(me) # growth rate is 0.3
    # print S
    print("add GSM constraint")
    # for i in metabolites:
    for i in S.keys():
        label = "mass_balance_%s"%i
        dot_S_v = pulp.lpSum([S[i][j] * v[j] for j in S[i].keys()])
        condition = dot_S_v == 0
        lp_prob += condition, label

    ###### cut in the genome
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
    print("add TU constraint")
    for gene,promoters in TUs.items():
        if gene in not_shared: continue
        gene = 'u_G_' + gene
        len_pro = len(promoters)
        print(gene, promoters)
        lp_prob += z[gene] - pulp.lpSum(z['u_G_'+j] for j in promoters) + \
                    (len_pro - 1) >= 0,'TU_all_'+gene
        for pro in promoters:
            pro = 'u_G_' + pro
            lp_prob += z[gene] - z[pro] <=0, 'TU_'+gene+'_'+pro

    ##### some overlapped region cannot be selected as the start of deletion
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
    GUROBI_CMD_OPTIONS = [('Threads', 8), ('TimeLimit', 1800), ('FeasibilityTol',1E-9),
                      ('OptimalityTol',1E-9),('IntFeasTol',1E-9),
                      ('MIPGapAbs', 0), ('MIPGap', 0), ('CliqueCuts', 2)]
    pulp_solver = pulp.solvers.GUROBI_CMD(path=None, keepFiles=0, mip=1, msg=0,
                            options=GUROBI_CMD_OPTIONS)
    # pulp_solver = pulp.solvers.CPLEX(path=None, keepFiles=0, mip=1,\
    #     msg=1, options=['mip tolerances mipgap 0', \
    #     'mip tolerances absmipgap 0', 'mip tolerances integrality 0',\
    #     'simplex tolerances optimality 1E-9',\
    #     'simplex tolerances feasibility 1E-9',], timelimit=1200)
    x_list = []
    y_list = []
    status = []
    def iterate_solve(lp_prob,iter_count):
        lp_prob.solve(pulp_solver)
        print("----------- " + str(iter_count) + " ------------")
        status.append(pulp.LpStatus[lp_prob.status])
        print("Status:", pulp.LpStatus[lp_prob.status])
        for v in lp_prob.variables():
            if "x_u_G_" in v.name and v.varValue == 1:               
                xname = v.name.replace("x_","")
                xname = xname.replace('_','-')
                xname = xname.replace("PM-","PM_")
                xname = xname.replace('u-','u_')
                xname = xname.replace('G-','G_')
                
                print(xname,v.name)
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

    for iter_count in range(1,2):
        lp_prob = iterate_solve(lp_prob,iter_count)
