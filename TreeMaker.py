# regular expression module
import re

class TreeTableNode:
    def __init__(self, name):
        self.data = ""
        self.name = str(name)
        self.children = []
        self.parent = []

    def add_child(self, obj):
        self.children.append(obj)

    def add_parent(self, obj):
        self.parent.append(obj)

    def get_parent(self):
        return self.parent[0]


class TreeTable:

    def __init__(self, inFile):
        fo = open(inFile, "r")
        fs = fo.read();
        processed = self.preprocess(fs)
        self.branches = self.branchFinder(processed)
        self.parents, self.children = self.parentChildFinder(self.branches)
        self.root
        # print for debugging
        print(processed)
        print(self.branches)
        print(self.parents)
        print(self.children)
        fo.close()

    def get_root(self):
        return self.root

    def get_tree(self):
        return parents, children

    # this removes branch length values
    # returns a Newick formatted string
    def preprocess(self, fileString):
        fileString = fileString.replace("\n","")
        sre = re.compile(':0\.\d+')
        fileString = sre.sub("",fileString)
        stringList = []
        pos = 0
        for char in fileString:
            if char in ['(', ')']:
                stringList.append(char)
                pos += 1
            elif char == ',':
                pos += 1
                continue
            else:
                name = ''
                while fileString[pos] not in ['(', ')', ',',';']:
                    name += fileString[pos]
                    pos += 1
                if name != '':
                    stringList.append(name)
        return stringList

    def branchFinder(self, processedString):
        branches = []
        opens = []
        bpos = 0
        
        node_counter = 1
        self.root = TreeTableNode(0)
        current_node = self.root

        for i, char in enumerate(processedString):
            if i == 0 or i == len(processedString) - 1:
                continue
            if char == '(':
                opens.append(i)
                node = TreeTableNode(node_counter)
                node.add_parent(current_node)
                current_node.add_child(node)
                current_node = node
                node_counter += 1
                
            elif char == ')':
                branches.append([])
                for letter in processedString[opens.pop()+1:i]:
                    if letter in [',', '(', ')']:
                        continue
                    else:
                        branches[bpos].append(letter)
                        if current_node != self.root:
                            current_node = current_node.get_parent()
                bpos += 1
            elif char == ',':
                continue
            else:
                branches.append([char])
                n = TreeTableNode(char)
                current_node.add_child(n)
                bpos += 1
        
        #print(root.children[1].children[0].name) # DEBUG
        return branches

    def parentChildFinder(self, branches):
        parents = []
        children = [[] for i in range(len(branches))]

        for i, branch in enumerate(branches):
            if i != len(branches) - 1:
                j = i + 1
                while not frozenset(branch) < frozenset(branches[j]):
                    if j == len(branches) - 1:
                        j = -1
                        break
                    j += 1
                parents.append(j)
                if j != -1:
                    children[j].append(i)
            else:
                parents.append(-1)
        return (parents, children)
