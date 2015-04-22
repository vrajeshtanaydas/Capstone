################################
#
#Treemaker converts a phylogenetic tree file into a tree data structure, 
#where each node in the tree contains:
#   data - used to store A T C or G, when a particular SNP is loaded      
#           into the tree using entropy calculator
#   name - stores the species name or branch number of the node depending  
#           on if the node is a leaf or not
#   children - stores the list of this nodes children
#   parent - stores this nodes parent
#
#The preprocess method removes unnecessary data from the tree file, and 
#creates a list containing "(", ")", and species names
#
#The branchfinder method reads through the processed list and builds the 
#initial tree
#
#the parentChildFinder method is used to generate the parents and 
#children for each branch
#
################################

# regular expression module
import re
import queue


# exception to demonstrate badly formatted tree (newick) file
class TreeFileError(Exception):
    pass


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

    def printTree(self, level=0):
        print ('\t' * level + repr(level) + self.data + '-' + repr(len(self.children)) + ' ' +self.name)
        for child in self.children:
            child.printTree(level+1)


class TreeTable:

    def __init__(self, inFile):
        self.fo = open(inFile, "r")
        fs = self.fo.read();

        # close the file when propogating the TreeFileError
        try:
            self.processed = self.preprocess(fs)
            self.branches = self.branchFinder(self.processed)
            self.parents, self.children = self.parentChildFinder(self.branches)
            self.fo.close()
        except TreeFileError as err:
            self.fo.close()
            raise(err)


    def get_root(self):
        return self.root

    def get_tree(self):
        return parents, children


    # preprocess -- removes branch length/distance values
    # -- removes branch distance values as well as commas 
    #    and any newlines
    # -- input: contents of a Newick file as a string
    # -- return: a list of parenthesis and genome names
    def preprocess(self, fileString):
        # strip any newlines from the file string
        fileString = fileString.replace("\n","")
        # remove any distances
        sre = re.compile(':0\.\d+')
        fileString = sre.sub("",fileString)
        # store parenthesis and genome names 
        processed_string = []
        # keeps track of current position in the processed file string
        pos = 0

        # iterate through each character in the processed file string
        for char in fileString:

            # we keep the parenthesis for determining tree structure
            if char in ['(', ')']:
                processed_string.append(char)
                pos += 1

            # don't include commas
            elif char == ',':
                pos += 1
                continue

            # store a genome name as a string
            else:
                name = ''
                try:
                    # continue adding chars to genome name if char isn't a
                    # delineator
                    while fileString[pos] not in ['(', ')', ',',';']:
                        name += fileString[pos]
                        pos += 1
                    # when a delineator is reached, append the genome name
                    if name != '':
                        processed_string.append(name)
                # occurs when an invalid character is last char of file string
                except IndexError:
                    self.fo.close()
                    raise TreeFileError("Third argument must be a properly" +
                                        " formatted Newick file" +
                                        "\n\t -- check last character in file")

        return processed_string


    # The branchfinder method reads through the processed list and builds the 
    # initial tree structure
    def branchFinder(self, processedString):
        branches = []
        # a list of open parens -- used for determining new tree levels
        opens = []
        # the current branch position
        bpos = 0
        # used for naming nodes that aren't leaf nodes
        node_identifier = 1
        # initialize the root node with name = 0
        self.root = TreeTableNode(0)
        # keeps track of which node currently at in tree structure 
        current_node = self.root

        # expected tokens are a comma, parens, or genome name
        for i, token in enumerate(processedString):

            # skip the first and last parenthesis
            # -- these don't determine the tree structure
            if i == 0 or i == len(processedString) - 1:
                
                # newick file must begin with open parens
                if i == 0 and token != '(':
                    raise TreeFileError("\nTree file is not properly formatted" +
                                        "\n\t-- must start with a \'(\'" +
                                        " character.\n")
                # newick file must end with close parens
                elif i == len(processedString) - 1 and token != ')':
                    raise TreeFileError("\nTree file is not properly formatted" +
                                        "\n\t-- must end with a \')\'" +
                                        " character\n")
                else:
                    continue

            # descending a tree level
            if token == '(':
                opens.append(i)
                # add a child node to the current node
                node = TreeTableNode(node_identifier)
                node.add_parent(current_node)
                current_node.add_child(node)
                # set the current node to the added child node
                # i.e. descend
                current_node = node
                node_identifier += 1
                
            # ascending a tree level(s)
            elif token == ')':
                # 
                branches.append([])
                try:
                    for token2 in processedString[opens.pop()+1:i]:
                        if token2 in [',', '(', ')']:
                            continue
                        else:
                            branches[bpos].append(token2)
                            if current_node != self.root:
                                current_node = current_node.parent[0]
                except IndexError:
                    raise TreeFileError("\nTree file is not properly formatted" +
                            "\n\t-- ensure number of open and close" +
                            " parentheses match\n")
                bpos += 1

            # add leaf nodes 
            else:
                branches.append([token])
                n = TreeTableNode(token)
                current_node.add_child(n)
                bpos += 1
        
        print(branches)
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
        print( (parents, children) )
        return (parents, children)

    
