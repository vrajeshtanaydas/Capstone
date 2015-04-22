#########################################################
#
#EntropyCalculator calculates the entropy value for a SNP or pair of SNPs
#by using the phylogenetic tree to calculate the probabilities of each value
#(A, T, C, G, AA, AT, ect.), using the methods demonstrated in:
#http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3049237/ as a reference.
#
#The probability of a given value is the number of time it occurs as a separate 
#group in the tree, divided by the total number of separate groups
#
#The entropy value is then:
#
#   sum( probability(value) * log(probability(value)) )   
#
#to use the entropy calculator, pass in a tree and the A, T , C, and G genome 
#sets from a SNP, (for two SNPs, pass in the second SNP's A, T, C, and G 
#groups as well)
#
#########################################################

import math
import TreeMaker

aGenomes = 0
tGenomes = 0
cGenomes = 0
gGenomes = 0

aGenomes2 = 0
tGenomes2 = 0
cGenomes2 = 0
gGenomes2 = 0

def visit(currentNode):
    childNum = len(currentNode.children)
    if childNum == 0:
        if currentNode.name in aGenomes:
            currentNode.data = 'A'
        elif currentNode.name in tGenomes:
            currentNode.data = 'T'
        elif currentNode.name in cGenomes:
            currentNode.data = 'C'
        elif currentNode.name in gGenomes:
            currentNode.data = 'G'
        else:
            return

        if currentNode.name in aGenomes2:
            currentNode.data += 'A'
        elif currentNode.name in tGenomes2:
            currentNode.data += 'T'
        elif currentNode.name in cGenomes2:
            currentNode.data += 'C'
        elif currentNode.name in gGenomes2:
            currentNode.data += 'G'

    else:
        for i in range(childNum):
            visit(currentNode.children[i])
        currentData = currentNode.children[0].data
        for j in range(childNum):
            if currentNode.children[j].data == 'X':
                continue
            elif currentNode.children[j].data != currentData \
                    and currentNode.children[j].data != 'X':
                currentNode.data = 'X'    
                break
            else:
                currentNode.data = currentNode.children[j].data
        

def incrementGroups(currentNode):
    global groups
    global a_group
    global c_group
    global t_group
    global g_group
    global aa_group
    global at_group
    global ac_group
    global ag_group
    global ta_group
    global tt_group
    global tc_group
    global tg_group
    global ca_group
    global ct_group
    global cc_group
    global cg_group
    global ga_group
    global gt_group
    global gc_group
    global gg_group
    groups+=1
    if currentNode.data == 'A':
        a_group += 1
    elif currentNode.data == 'T':
        t_group += 1
    elif currentNode.data == 'C':
        c_group += 1
    elif currentNode.data == 'G':
        g_group += 1
    elif currentNode.data == 'AA':
        aa_group += 1
    elif currentNode.data == 'AT':
        at_group += 1
    elif currentNode.data == 'AC':
        ac_group += 1
    elif currentNode.data == 'AG':
        ag_group += 1
    elif currentNode.data == 'TA':
        ta_group += 1
    elif currentNode.data == 'TT':
        tt_group += 1
    elif currentNode.data == 'TC':
        tc_group += 1
    elif currentNode.data == 'TG':
        tg_group += 1
    elif currentNode.data == 'CA':
        ca_group += 1
    elif currentNode.data == 'CT':
        ct_group += 1
    elif currentNode.data == 'CC':
        cc_group += 1
    elif currentNode.data == 'CG':
        cg_group += 1
    elif currentNode.data == 'GA':
        ga_group += 1
    elif currentNode.data == 'GT':
        gt_group += 1
    elif currentNode.data == 'GC':
        gc_group += 1
    elif currentNode.data == 'GG':
        gg_group += 1


#Calculate Value determines the entropy value
def CalculateValue(treeTable, currentNode):
    global groups
    global a_group
    global c_group
    global t_group
    global g_group
    global aa_group
    global at_group
    global ac_group
    global ag_group
    global ta_group
    global tt_group
    global tc_group
    global tg_group
    global ca_group
    global ct_group
    global cc_group
    global cg_group
    global ga_group
    global gt_group
    global gc_group
    global gg_group

    childCharacter = []
    childNum = len(currentNode.children)
    if childNum == 0:
        pass
    elif currentNode.data == 'X' and currentNode != treeTable.get_root():
        childCharacter.append(currentNode.parent[0].data)
    elif currentNode == treeTable.get_root() and currentNode.data != 'X':
        childCharacter.append(currentNode.data)
        incrementGroups(currentNode)
    else:
        childCharacter.append(currentNode.data)
    if childNum != 0:
        for i in range(childNum):
            if currentNode.children[i].data in childCharacter or \
                    currentNode.children[i].data == 'X':
                pass
            else:
                childCharacter.append(currentNode.children[i].data)
                incrementGroups(currentNode.children[i])
            CalculateValue(treeTable,currentNode.children[i])


def CalculateEntropy():
    global groups
    global a_group
    global t_group
    global c_group
    global g_group
    global aa_group
    global at_group
    global ac_group
    global ag_group
    global ta_group
    global tt_group
    global tc_group
    global tg_group
    global ca_group
    global ct_group
    global cc_group
    global cg_group
    global ga_group
    global gt_group
    global gc_group
    global gg_group
    p_a = 0
    p_t = 0
    p_c = 0
    p_g = 0
    p_aa = 0
    p_at = 0
    p_ac = 0
    p_ag = 0
    p_ta = 0
    p_tt = 0
    p_tc = 0
    p_tg = 0
    p_ca = 0
    p_ct = 0
    p_cc = 0
    p_cg = 0
    p_ga = 0
    p_gt = 0
    p_gc = 0
    p_gg = 0
    if a_group != 0:
        p_a = a_group / groups
        p_a *= math.log2(p_a)
    if t_group != 0:
        p_t = t_group / groups
        p_t *= math.log2(p_t)
    if c_group != 0:
        p_c = c_group / groups
        p_c *= math.log2(p_c)
    if g_group != 0:
        p_g = g_group / groups
        p_g *= math.log2(p_g)
    if aa_group != 0:
        p_aa = aa_group / groups
        p_aa *= math.log2(p_aa)
    if at_group != 0:
        p_at = at_group / groups
        p_at *= math.log2(p_at)
    if ac_group != 0:
        p_ac = ac_group / groups
        p_ac *= math.log2(p_ac)
    if ag_group != 0:
        p_ag = ag_group / groups
        p_ag *= math.log2(p_ag)
    if ta_group != 0:
        p_ta = ta_group / groups
        p_ta *= math.log2(p_ta)
    if tt_group != 0:
        p_tt = tt_group / groups
        p_tt *= math.log2(p_tt)
    if tc_group != 0:
        p_tc = tc_group / groups
        p_tc *= math.log2(p_tc)
    if tg_group != 0:
        p_tg = tg_group / groups
        p_tg *= math.log2(p_tg)
    if ca_group != 0:
        p_ca = ca_group / groups
        p_ca *= math.log2(p_ca)
    if ct_group != 0:
        p_ct = ct_group / groups
        p_ct *= math.log2(p_ct)
    if cc_group != 0:
        p_cc = cc_group / groups
        p_cc *= math.log2(p_cc)
    if cg_group != 0:
        p_cg = cg_group / groups
        p_cg *= math.log2(p_cg)
    if ga_group != 0:
        p_ga = ga_group / groups
        p_ga *= math.log2(p_ga)
    if gt_group != 0:
        p_gt = gt_group / groups
        p_gt *= math.log2(p_gt)
    if gc_group != 0:
        p_gc = gc_group / groups
        p_gc *= math.log2(p_gc)
    if gg_group != 0:
        p_gg = gg_group / groups
        p_gg *= math.log2(p_gg)

    entropy = -1 * (p_a + p_t + p_c + p_g +
            p_aa + p_at + p_ac + p_ag +
            p_ta + p_tt + p_tc + p_tg +
            p_ca + p_ct + p_cc + p_cg +
            p_ga + p_gt + p_gc + p_gg)
    return entropy


def main(treeTable, aSet, tSet, cSet, gSet, aSet2=[], tSet2=[],cSet2=[],gSet2=[]):
    global aGenomes 
    aGenomes = aSet
    global tGenomes
    tGenomes = tSet
    global cGenomes 
    cGenomes = cSet 
    global gGenomes 
    gGenomes = gSet

    global aGenomes2
    aGenomes2 = aSet2
    global tGenomes2
    tGenomes2 = tSet2
    global cGenomes2
    cGenomes2 = cSet2
    global gGenomes2
    gGenomes2 = gSet2



    visit(treeTable.get_root())
    global groups
    global a_group
    global t_group
    global c_group
    global g_group
    global aa_group
    global at_group
    global ac_group
    global ag_group
    global ta_group
    global tt_group
    global tc_group
    global tg_group
    global ca_group
    global ct_group
    global cc_group
    global cg_group
    global ga_group
    global gt_group
    global gc_group
    global gg_group
    groups = 0
    a_group = 0
    t_group = 0
    c_group = 0
    g_group = 0
    aa_group = 0
    at_group = 0
    ac_group = 0
    ag_group = 0
    ta_group = 0
    tt_group = 0
    tc_group = 0
    tg_group = 0
    ca_group = 0
    ct_group = 0
    cc_group = 0
    cg_group = 0
    ga_group = 0
    gt_group = 0
    gc_group = 0
    gg_group = 0

    CalculateValue(treeTable, treeTable.get_root())

    entropy = CalculateEntropy()
    return entropy
