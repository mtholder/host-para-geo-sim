#!/usr/bin/env python
import sys
import dendropy
import os
if len(sys.argv) < 2:
    sys.exit('''Error: expected a filepath as the first (and only argument).
    
Program:
    summary_stats.py  Copyright (C) 2014 Mark T. Holder mtholder@gmail.com
    This program comes with ABSOLUTELY NO WARRANTY. See LICENSE.txt for details.

    Calculator of some simple summary statistics that quantify the degree of 
    host-parasite co-phylogeny. Parses the quirky output of sim-host-parasite.py
    as its input.

Prerequisites:
    DendroPy

Usage:
    python summary_stats.py </path/to/output/of/sim-host-parasite.py>

''')
PRUNING_SINGLE_OUTGROUP = os.environ.get('PRUNING_SINGLE_OUTGROUP', '0') == '1' # should make this more flexible for >1 outgroup...
REVERSE_MAPPING = os.environ.get('REVERSE_MAPPING', '0') == '1' # should make this more flexible for >1 outgroup...
sys.stderr.write('PRUNING_SINGLE_OUTGROUP = {}, REVERSE_MAPPING={}\n'.format(PRUNING_SINGLE_OUTGROUP, REVERSE_MAPPING))
inpfn = sys.argv[1]
inp = open(inpfn)
lines = inp.readlines()
h_newick = lines[1]
p_newick = lines[2]
p2h = {}
for line in lines[3:]:
    ls = line.strip().split()
    p = ls[0]
    hl = ls[1:]
    p2h[p] = hl

h_taxa = dendropy.TaxonSet()
h_tree = dendropy.Tree.get_from_string(h_newick, schema='newick', taxon_set=h_taxa)
p_taxa = dendropy.TaxonSet()
p_tree = dendropy.Tree.get_from_string(p_newick, schema='newick', taxon_set=p_taxa)
if PRUNING_SINGLE_OUTGROUP:
    for tree_taxa in [(h_tree, h_taxa), (p_tree, p_taxa)]:
        _tree, _taxa = tree_taxa
        print _tree
        _r = _tree.seed_node
        _c = _r.child_nodes()
        assert len(_c) == 2
        _ingroup, _outgroup = _c
        if _c[1].child_nodes():
            assert not _c[0].child_nodes()
            _outgroup, _ingroup = _c
        else:
            assert _c[0].child_nodes()
        _outgroup = _outgroup.taxon
        _tree.seed_node = _ingroup
        _ingroup.edge.tail_node = None
        _ingroup.edge.length = 0.0
        _ingroup.edge.rootedge = True
        _oind = None
        for n, t in enumerate(_taxa):
            if t is _outgroup:
                _oind = n
                break
        assert _oind is not None
        del _taxa[_oind]
        print _tree
if REVERSE_MAPPING:
    h2p = {}
    for h in h_taxa:
        h2p[h.label] = []
    h_tree, p_tree = p_tree, h_tree
    h_taxa, p_taxa = p_taxa, h_taxa
    for k, v in p2h.items():
        for h in v:
            if h in h2p:
                h2p[h].append(k)
    p2h = h2p
    del h2p
h_lists_list = []
taxon_p2h = {}
for k, v in p2h.items():
    val_list = [h_taxa.get_taxon(label=i) for i in v if h_taxa.get_taxon(label=i) is not None]
    h_lists_list.append(val_list)
    taxon_p2h[p_taxa.get_taxon(label=k)] = val_list

def calculate_mean_sum_of_tree_length_covered(h_lists_list, h_tree):
    t0 = 0.0
    c = 0
    for v in h_lists_list:
        c += 1
        if len(v) > 1:
            mrca_root = h_tree.mrca(taxa=v)
            d = {}
            for n in v:
                curr_nd = h_tree.find_node_with_taxon(lambda x: x == n)
                while (curr_nd is not mrca_root) and (curr_nd not in d):
                    d[curr_nd] = curr_nd.edge.length
                    curr_nd = curr_nd.parent_node
            s = sum([v for v in d.values()])
            t0 += s
    return t0/c

h_len = h_tree.length()
p_tree.calc_node_ages(check_prec=0.00001)
p_age = p_tree.seed_node.age
tl_stats = []
t0 = calculate_mean_sum_of_tree_length_covered(h_lists_list, h_tree)
tl_stats.append((0.0, t0))
for nd in p_tree.postorder_node_iter():
    if nd.is_leaf():
        htl = taxon_p2h[nd.taxon]
        nd.hts = set(htl)
    else:
        nd.hts = set()
        for c in nd.child_nodes():
            nd.hts.update(c.hts)

n_quantiles = 100
for i in range(1, n_quantiles):
    proportion = i/float(n_quantiles)
    crit = proportion*p_age
    h_lists_list = []
    for e in p_tree.postorder_edge_iter():
        if e.head_node is not p_tree.seed_node:
            if (e.tail_node.age >= crit) and (e.head_node.age < crit):
                h_lists_list.append(e.head_node.hts)
    tp = calculate_mean_sum_of_tree_length_covered(h_lists_list, h_tree)/h_len
    tl_stats.append((proportion, tp))

h_lists_list=[p_tree.seed_node.hts]
tmax = calculate_mean_sum_of_tree_length_covered(h_lists_list, h_tree)/h_len
tl_stats.append((1.0, tmax))


tr_stats = []
for p, c in enumerate(p_tree.seed_node.child_nodes()):
    h_lists_list=[c.hts]
    tp = calculate_mean_sum_of_tree_length_covered(h_lists_list, h_tree)/h_len
    tr_stats.append((p, tp))

for p, tp in tl_stats:
    print 't{p:3.2f} = {s:f}'.format(p=p, s=tp)
for p, tp in tr_stats:
    print 'tr{p:d} = {s:f}'.format(p=p, s=tp)
for p, tp in tl_stats:
    print 'r{p:3.2f} = {s:f}'.format(p=p, s=tp/tmax)
for p, tp in tr_stats:
    print 'rr{p:d} = {s:f}'.format(p=p, s=tp/tmax)



# Calculator of some simple summary statistics that quantify the degree of 
#   host-parasite co-phylogeny. Parses the quirky output of sim-host-parasite.py
#   as its input.
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
