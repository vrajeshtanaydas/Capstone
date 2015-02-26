import math
import TreeTableMaker

aGenomes = 0
tGenomes = 0
cGenomes = 0
gGenomes = 0

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
            print(currentNode.name)
        return
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
        
#Calculate Value determines the entropy value
def CalculateValue(treeTable, currentNode):
    childCharacter = []
    childNum = len(currentNode.children)
    #print("\n\tCurrent: " + str(currentNode.data) + "\t NumChildren: " + \
    #        str(childNum))
    if childNum == 0:
        pass
    elif currentNode.data == 'X' and currentNode != treeTable.get_root():
        childCharacter.append(currentNode.get_parent().data)
    else:
        childCharacter.append(currentNode.data)
    if childNum != 0:
        for i in range(childNum):
            if currentNode.children[i].data in childCharacter or \
                    currentNode.children[i].data == 'X':
                pass
            else:
                childCharacter.append(currentNode.children[i].data)
                global groups
                groups+=1
                if currentNode.children[i].data == 'A':
                    global a_group
                    a_group += 1
                elif currentNode.children[i].data == 'T':
                    global t_group
                    t_group += 1
                elif currentNode.children[i].data == 'C':
                    global c_group
                    c_group += 1
                elif currentNode.children[i].data == 'G':
                    global g_group
                    g_group += 1

            CalculateValue(treeTable,currentNode.children[i])

def BFS(currentNode):
    print("\t" + currentNode.data)
    childNum = len(currentNode.children)
    for i in range(childNum):
        BFS(currentNode.children[i])


def main(treeTable, aSet, tSet, cSet, gSet):
    global aGenomes 
    aGenomes = aSet
    global tGenomes
    tGenomes = tSet
    global cGenomes 
    cGenomes = cSet 
    global gGenomes 
    gGenomes = gSet
    visit(treeTable.get_root())
    #treeTable.get_root().data = 'X'
    print("NEW TABLE")
    BFS(treeTable.get_root())
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
    print(groups)
    print(a_group)
    print(t_group)
    print(c_group)
    print(g_group)
    return treeTable
