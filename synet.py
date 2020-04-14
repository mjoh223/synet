from collections import defaultdict

IDs = ['a','b', 'c','d','e']
loci_di = default_dict('key':'values')
print(loci_di)
class Locus(object):
    def __init__(self, anchor):
        self.anchor = anchor

for anchor in IDs:
    a = Locus(anchor)
    print(a.anchor)
