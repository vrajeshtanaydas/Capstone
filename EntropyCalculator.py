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
    groups+=1
    if currentNode.data == 'A':
        a_group += 1
    elif currentNode.data == 'T':
        t_group += 1
    elif currentNode.data == 'C':
        c_group += 1
    elif currentNode.data == 'G':
        g_group += 1

#Calculate Value determines the entropy value
def CalculateValue(treeTable, currentNode):
    global groups
    global a_group
    global c_group
    global t_group
    global g_group

    childCharacter = []
    childNum = len(currentNode.children)
    #print("\n\tCurrent: " + str(currentNode.data) + "\t NumChildren: " + \
    #        str(childNum))
    if childNum == 0:
        pass
    elif currentNode.data == 'X' and currentNode != treeTable.get_root():
        childCharacter.append(currentNode.parent[0].data)
    elif currentNode == treeTable.get_root() and currentNode.data != 'X':
       # print(currentNode.data)
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
               # print(currentNode.data)
                incrementGroups(currentNode.children[i])

            CalculateValue(treeTable,currentNode.children[i])


def CalculateEntropy():
    global groups
    global a_group
    global t_group
    global c_group
    global g_group
    p_a = 0
    p_t = 0
    p_c = 0
    p_g = 0
    if a_group != 0:
        p_a = int(a_group) / groups
        p_a *= math.log10(p_a)
    if t_group != 0:
        p_t = int(t_group) / groups
        p_t *= math.log10(p_t)
    if c_group != 0:
        p_c = int(c_group) / groups
        p_c *= math.log10(p_c)
    if g_group != 0:
        p_g = int(g_group) / groups
        p_g *= math.log10(p_g)

    entropy = -1 * (p_a + p_t + p_c + p_g)
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
    #treeTable.get_root().data = 'X'
   # print("NEW TABLE")
    global groups
    global a_group
    global t_group
    global c_group
    global g_group
    groups = 0
    a_group = 0
    t_group = 0
    c_group = 0
    g_group = 0
    CalculateValue(treeTable, treeTable.get_root())

    entropy = CalculateEntropy()
    return entropy
