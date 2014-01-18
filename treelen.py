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
        

# Simple text based tool for summing branch lengths. This is very dumb and would
#   be confused by names with : in them. It uses regex to sum the numbers after
#   each colon.
#   Copyright (C) 2014 Mark T. Holder mtholder@gmail.com
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
