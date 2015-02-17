# regular expression module
import re

class TreeTable(object):

    def __init__(self, inFile):
        #fo = open("testTree.txt", "r")
        fo = open(inFile, "r")
        fs = fo.read();
        processed = self.preprocess(fs)
        self.branches = self.branchFinder(processed)
        self.parents, self.children = self.parentChildFinder(self.branches)
        # print for debugging
        print(processed)
        print(self.branches)
        print(self.parents)
        print(self.children)
        fo.close()


    def get_siblings(self, species):
        pass

    def get_tree(self):
        return self

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

        for i, char in enumerate(processedString):
            if i == 0 or i == len(processedString) - 1:
                continue
            if char == '(':
                opens.append(i)
            elif char == ')':
                branches.append([])
                for letter in processedString[opens.pop()+1:i]:
                    if letter in [',', '(', ')']:
                        continue
                    else:
                        branches[bpos].append(letter)
                bpos += 1
            elif char == ',':
                continue
            else:
                branches.append([char])
                bpos += 1
        return branches

    def parentChildFinder(self, branches):
        parents = []
        children = [[] for i in range(len(branches))]

        for i, branch in enumerate(branches):
            if i != len(branches) - 1:
                j = i + 1
                while not set(branch) < set(branches[j]):
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
