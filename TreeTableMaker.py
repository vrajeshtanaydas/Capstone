# regular expression module
import re

# this removes branch length values
# returns a Newick formatted string
def preprocess(fileString):
    fileString = fileString.replace("\n","")
    sre = re.compile(':0\.\d+')
    return sre.sub("",fileString)

fo = open("test.txt", "r")

fs = fo.read();
processed = preprocess(fs)

print(processed)

fo.close()


  
    










