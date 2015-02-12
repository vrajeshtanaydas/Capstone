# regular expression module
import re

# this removes branch length values
# returns a Newick formatted string
def preprocess(fileString):
    fileString = fileString.replace("\n","")
    sre = re.compile(':0\.\d+')
    fileString = sre.sub(" ",fileString)
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
            while fileString[pos] not in ['(', ')', ',']:
                name += fileString[pos]
                pos += 1
            if name != '':
                stringList.append(name)
    return fileString

def branchFinder(fileString):
    branches = []
    opens = []
    bpos = 0

    for i, char in enumerate(fileString):
        if i == 0 or i == len(fileString) - 1:
            continue
        if char == '(':
            opens.append(i)
        elif char == ')':
            branches.append([])
            for letter in fileString[opens.pop()+1:i]:
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

def parentChildFinder(branches):
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

fo = open("test.txt", "r")

fs = fo.read();
processed = preprocess(fs)
branches = branchFinder(processed)
parents, children = parentChildFinder(branches)

print(processed)
print(branches)
print(parents)
print(children)

fo.close()


  
    










