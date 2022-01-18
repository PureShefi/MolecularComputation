from mml import *

lab = Lab()
lab.samples.append(DNA("aaa", "ttt"))
lab.samples.append(DNA("aaa", "ttt"))
lab.samples.append(DNA("aaa", "ttt"))
lab.samples.append(DNA("ggg", "ccc"))

print(lab.samples)

lab.amplify(["a", "t"], 2)
print(lab.samples)

lab.samples.append(DNA("ggg", ""))
lab.samples.append(DNA("    c", "cccg"))
lab.merge()

print(lab.samples)