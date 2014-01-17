#!/usr/bin/env python
import sys
import dendropy

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
h_lists_list = []
taxon_p2h = {}
for k, v in p2h.items():
    val_list = [h_taxa.get_taxon(label=i) for i in v]
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
                while (curr_nd is not mrca_root) and (curr_nd.taxon not in d):
                    d[curr_nd.taxon] = curr_nd.edge.length
                    curr_nd = curr_nd.parent_node
            s = sum([v for v in d.values()])
            t0 += s
    return t0/c

h_len = h_tree.length()
p_tree.calc_node_ages()
p_age = p_tree.seed_node.age
t0 = calculate_mean_sum_of_tree_length_covered(h_lists_list, h_tree)
print 't0.00 =', t0/h_len
for nd in p_tree.postorder_node_iter():
    if nd.is_leaf():
        htl = taxon_p2h[nd.taxon]
        nd.hts = set(htl)
    else:
        nd.hts = set()
        for c in nd.child_nodes():
            nd.hts.update(c.hts)

for i in range(1,10):
    proportion = i/10.0
    crit = proportion*p_age
    h_lists_list = []
    for e in p_tree.postorder_edge_iter():
        if e.head_node is not p_tree.seed_node:
            if (e.tail_node.age >= crit) and (e.head_node.age < crit):
                h_lists_list.append(e.head_node.hts)
    tp = calculate_mean_sum_of_tree_length_covered(h_lists_list, h_tree)/h_len
    print 't{p:3.2f} = {s:f}'.format(p=proportion, s=tp)
