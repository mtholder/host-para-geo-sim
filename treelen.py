#!/usr/bin/env python
import sys, re
pat = re.compile(r":([-.0-9Ee]+)")
try:
    tree = sys.argv[1]
    if tree.strip().startswith('('):
        m = pat.findall(tree)
        print sum([float(i) for i in m])
    else:
        f = open(tree, 'rU')
        for line in f:
            m = pat.findall(line)
            print sum([float(i) for i in m])
except:
    for line in sys.stdin:
        m = pat.findall(line)
        print sum([float(i) for i in m])
        
