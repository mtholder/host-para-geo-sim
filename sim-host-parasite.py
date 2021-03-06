#!/usr/bin/env python

_EMPTY_SET = frozenset()
SILENT_MODE = False
def debug(m):
    global DEBUGGING_OUTPUT
    if DEBUGGING_OUTPUT:
        sys.stderr.write('debug: ')
        sys.stderr.write(m)
        sys.stderr.write('\n')

def warn(m):
    global SILENT_MODE
    if not SILENT_MODE:
        sys.stderr.write('warn: ')
        sys.stderr.write(m)
        sys.stderr.write('\n')
def error(m):
    sys.stderr.write('Error: ')
    sys.stderr.write(m)
    sys.stderr.write('\n')

class Locality(object):
    def __init__(self, latitude, longitude):
        self.lat, self.lon = latitude, longitude
        self._h2p = {}
        self._p2h = {}
        self.h_prev_time = {}
        # north, east, south, and west
        self._n_neighbor = None
        self._e_neighbor = None
        self._s_neighbor = None
        self._w_neighbor = None
        self._hash = hash((latitude, longitude))
    def __hash__(self):
        return self._hash
    def _remove_h_from_p2h(self, h):
        op = self._h2p.get(h)
        if op:
            for p in op:
                ps = self._p2h[p]
                ps.remove(h)
                if not ps:
                    del self._p2h[p]
                    p.geo_range.remove(self)
    def _remove_p_from_h2p(self, p):
        oh = self._p2h.get(p)
        if oh:
            for h in oh:
                ps = self._h2p[h]
                ps.remove(p)
                #if not ps:
                #    del self._h2p[h]
                #    h.geo_range.remove(self)
    def _add_h_to_p2h(self, h, p_set):
        for p in p_set:
            ps = self._p2h.get(p)
            if ps is None:
                self._p2h[p] = set([h])
            else:
                ps.add(h)
            p.geo_range.add(self)
    def set_para_range(self, h, p_set):
        self._remove_h_from_p2h(h)
        self._h2p[h] = p_set
        self._add_h_to_p2h(h, p_set)
    def remove_host(self, h):
        self._remove_h_from_p2h(h)
        del self._h2p[h]
    def remove_parasite_broken_range(self, p):
        self._remove_p_from_h2p(p)
        if p in self._p2h:
            del self._p2h[p]
    def remove_para_from_host(self, para, h):
        p_set = self._h2p[h]
        p_set.remove(para)
        #if not p_set:
        #    del self._h2p[h]
        #    h.geo_range.remove(self)
        h_set = self._p2h[para]
        h_set.remove(h)
        if not h_set:
            del self._p2h[para]
            para.geo_range.remove(self)
    def remove_host_if_present(self, h):
        if h in self._h2p:
            self.remove_host(h)
            h.geo_range.remove(self)
    def remove_para_from_host_if_present(self, para, h):
        hs = self._p2h.get(para)
        if hs and h in hs:
            self.remove_para_from_host(para, h)
    def get_para_for_host(self, h):
        global _EMPTY_SET
        return  self._h2p.get(h, _EMPTY_SET)
    def get_hosts_for_para(self, para):
        global _EMPTY_SET
        return  self._p2h.get(para, _EMPTY_SET)
    def add_para_for_h_host_range_broken(self, h, p_set):
        old_p_set = self._h2p.setdefault(h, set())
        #assert(self in h.geo_range)
        for p in p_set:
            if p not in old_p_set:
                old_p_set.add(p)
                self._p2h.setdefault(p, set()).add(h)
                p.geo_range.add(self)
    def replace_host(self, old_host, new_host):
        '''transfer all of the parasites from old_host to new_host (removing old_host).'''
        ov = self.get_para_for_host(old_host)
        if ov:
            self._h2p[new_host] = ov
            del self._h2p[old_host]
            for p in ov:
                hs = self._p2h[p]
                hs.remove(old_host)
                hs.add(new_host)
    def replace_parasite(self, old_para, new_para):
        h_set = self._p2h.get(old_para)
        if h_set is None:
            return
        del self._p2h[old_para]
        self._p2h.setdefault(new_para, set()).update(h_set)
        for h in h_set:
            hps = self._h2p[h]
            hps.remove(old_para)
            hps.add(new_para)
        old_para.geo_range.remove(self)
        new_para.geo_range.add(self)

    def set_neighbors(self, n):
        self._n_neighbor, self._e_neighbor, self._s_neighbor, self._w_neighbor = n
    def get_neighbors(self):
        return (self._n_neighbor, self._e_neighbor, self._s_neighbor, self._w_neighbor)
    neighbors = property(get_neighbors, set_neighbors)
    def get_east(self):
        return self._e_neighbor
    east = property(get_east)
    def get_west(self):
        return self._w_neighbor
    west = property(get_west)
    def get_north(self):
        return self._n_neighbor
    north = property(get_north)
    def get_south(self):
        return self._s_neighbor
    south = property(get_south)
    def __str__(self):
        return 'Locality({x:d}, {y:d})'.format(x=self.lon, y=self.lat)

class Grid(object):
    is_tube = False
    def __init__(self, width, height):
        assert(width > 0)
        assert(isinstance(width, int))
        assert(height > 0)
        assert(isinstance(height, int))
        self._w = width
        self._h = height
        self.rows = []
        self.cols = []
        for lat in range(height):
            r = [Locality(latitude=lat, longitude=i) for i in range(width)]
            self.rows.append(r)
        for j in range(width):
            c = [self.rows[i][j] for i in range(height)]
            self.cols.append(c)
        for i in range(height):
            for j in range(width):
                north = self.rows[i - 1][j] if i > 0 else None
                south = self.rows[i + 1][j] if i < (height -1) else None
                if j < (width - 1):
                    east = self.rows[i][j + 1]
                elif Grid.is_tube:
                    east = self.rows[i][0]
                else:
                    east = None
                if j > 0:
                    west = self.rows[i][j - 1]
                elif Grid.is_tube:
                    west = self.rows[i][width - 1]
                else:
                    west = None
                locality = self.rows[i][j]
                locality.neighbors = (north, east, south, west)
    def get(self, x, y):
        return self.cols[x % self._w][y]
    def gen_rand_range(self, max_ns_extent, max_ew_extent, inc_prob, rng):
        c_x = rng.randrange(self.width)
        c_y = rng.randrange(self.height)
        n_d = rng.randrange(max_ns_extent // 2)
        s_d = rng.randrange(max_ns_extent // 2)
        e_d = rng.randrange(max_ew_extent // 2)
        w_d = rng.randrange(max_ns_extent // 2)
        if Grid.is_tube:
            max_x = c_x + e_d
            min_x = c_x - w_d
        else:
            max_x = min(c_x + e_d, self.width - 1)
            min_x = max(c_x - w_d, 0)
        max_y = min(c_y + n_d, self.height - 1)
        min_y = max(c_y - s_d, 0)
        end_x = max_x + 1
        end_y = max_y + 1
        r = set()
        for x in range(min_x, end_x):
            c = self.get(x, c_y)
            r.add(c)
        for y in range(min_y, end_y):
            c = self.get(c_x, y)
            r.add(c)
        # SE quadrant
        for x in range(c_x + 1, end_x):
            for y in range(c_y + 1, end_y):
                c = self.get(x, y)
                if (c.north in r) or (c.west in r):
                    if h_rng.random() < inc_prob:
                        r.add(c)
        # NE quadrant
        for x in range(c_x + 1, end_x):
            for y in range(c_y, min_y -1, -1):
                c = self.get(x, y)
                if (c.south in r) or (c.west in r):
                    if h_rng.random() < inc_prob:
                        r.add(c)
        # SW quadrant
        for x in range(c_x - 1, min_x - 1, -1):
            for y in range(c_y + 1, end_y):
                c = self.get(x, y)
                if (c.north in r) or (c.east in r):
                    if h_rng.random() < inc_prob:
                        r.add(c)
        # NE quadrant
        for x in range(c_x - 1, min_x - 1, -1):
            for y in range(c_y, min_y-1, -1):
                c = self.get(x, y)
                if (c.south in r) or (c.east in r):
                    if h_rng.random() < inc_prob:
                        r.add(c)
        return r
    def get_width(self):
        return self._w
    width = property(get_width)
    def get_height(self):
        return self._h
    height = property(get_height)
    def write_range(self, o, geo_range):
        m = 2
        presence = unichr(0x2588)*m
        absence = ' '*m
        side_border = unichr(0x2551)
        horiz_border = unichr(0x2550)
        o.write(unichr(0x2554) + (self.width*m)*horiz_border + unichr(0x2557) + '\n')
        for r in self.rows:
            o.write(side_border)
            for j in range(self.width):
                c = r[j]
                if c in geo_range:
                    o.write(presence)
                else:
                    o.write(absence)
            o.write(side_border + '\n')
        o.write(unichr(0x255A) + (self.width*m)*horiz_border + unichr(0x255D) + '\n')

def find_periphery(r):
    p = set()
    for loc in r:
        n = 0
        ns = 0
        if loc.north not in r:
            ns += 1
            n += 1
        if loc.east not in r:
            n += 1
        if loc.south not in r:
            ns += 1
            n += 1
        if loc.west not in r:
            n += 1
        if n < 2:
            p.add(loc)
        elif n == 2:
            if ns != 0 and ns != 2:
                p.add(loc)
    return p

def find_surrounding(r):
    s = {}
    for loc in r:
        for n in loc.neighbors:
            if (n is not None) and (n not in r):
                s.setdefault(n, []).append(loc)
    return s

_HOST_LINEAGE_COUNTER = 0
_PARASITE_LINEAGE_COUNTER = 0

class BaseLineage(object):
    def __init__(self, start_time, parent):
        self.birth_time = start_time
        self.death_time = None
        self.parent = parent
        if parent:
            parent.children.append(self)
        self.children = []
        self._on_extant_trace = False
    def __hash__(self):
        return hash(self.index)
    def _get_des_extant_fork(self):
        ec = [self]
        while len(ec) == 1:
            nd = ec[0]
            ec = [i for i in nd.children if i._on_extant_trace]
        return nd, ec

class HostLineage(BaseLineage):
    def __init__(self, start_time, geo_range, parent=None):
        global _HOST_LINEAGE_COUNTER, CURR_TIME
        self.index = _HOST_LINEAGE_COUNTER
        _HOST_LINEAGE_COUNTER += 1
        BaseLineage.__init__(self, start_time, parent)
        self.geo_range = set(geo_range) # make a shallow copy of the localities
        for loc in self.geo_range:
            loc.set_para_range(self, set())
        if parent:
            self._mrca_times = dict(parent._mrca_times)
            for k, v in parent._mrca_times.items():
                k._mrca_times[self] = v
            self._mrca_times[parent] = CURR_TIME
            parent._mrca_times[self] = CURR_TIME
        else:
            self._mrca_times = {}
    def __str__(self):
        return 'HostLineage.index={i:d}'.format(i=self.index)

    def disperse(self, rng):
        global GRID, H_LOCATION_EXTINCTION_PROB, H_RANGE_EXPANSION_PROB
        p = find_periphery(self.geo_range)
        to_del = [loc for loc in p if rng.random() < H_LOCATION_EXTINCTION_PROB]
        for loc in to_del:
            loc.remove_host_if_present(self)
            SANITY_CHECK()
        p = find_surrounding(self.geo_range)
        added = []
        for k, v in p.items():
            if rng.random() < H_RANGE_EXPANSION_PROB:
                src = rng.sample(v, 1)[0]
                added.append((k, src))
        if added:
            for loc, src in added:
                para = src.get_para_for_host(self)
                if para is not None:
                    loc.add_para_for_h_host_range_broken(self, para)
            self.geo_range.update(set([i[0] for i in added]))
        SANITY_CHECK()
    def collect_para2range(self):
        d = {}
        for loc in self.geo_range:
            for p in loc.get_para_for_host(self):
                pr = d.setdefault(p, set())
                pr.add(loc)
        return d
    def close_enough(self, h_set):
        global host_infection_patristic_threshold, CURR_TIME
        if PHYLOGENETIC_HOST_JUMPS:
            for h in h_set:
                mrca_threshold = CURR_TIME - host_infection_patristic_threshold/2
                if self._mrca_times[h] >= mrca_threshold:
                    return True
            return False
        else:
            return True
    def debug_write(self, out):
        out.write('HostLineage')
        pi = str(self.parent.index) if self.parent is not None else 'None'
        out.write(' index = {i:d}\n parent.index = {p}\n'.format(i=self.index, p=pi))
        out.write(' geo_range:\n')
        for r in self.geo_range:
            p = [i.index for i in r.get_para_for_host(self)]
            p.sort()
            out.write('    {loc}: parasites = {pr}\n'.format(loc=str(r), pr=repr(p)))
        ot = [(v, k.index) for k, v in self._mrca_times.items()]
        ot.sort()
        out.write(' time to MRCAs:\n')
        for t, i in ot:
            out.write('    index = {i:d}: time = {t:d}\n'.format(i=i, t=t))

class ParasiteLineage(BaseLineage):
    def __init__(self, start_time, geo_range, host=None, parent=None):
        global _PARASITE_LINEAGE_COUNTER
        self.index = _PARASITE_LINEAGE_COUNTER
        _PARASITE_LINEAGE_COUNTER += 1
        BaseLineage.__init__(self, start_time, parent)
        self.geo_range = set()
        if host is not None:
            ls = [self]
            for loc in geo_range:
                loc.add_para_for_h_host_range_broken(host, ls)
    def __str__(self):
        return 'ParasiteLineage.index={i:d}'.format(i=self.index)
    def get_full_range(self):
        return self.geo_range
    def get_hosts(self):
        hs = set()
        for r in self.geo_range:
            hs.update(r.get_hosts_for_para(self))
        return hs
    def disperse(self, rng, tree):
        global GRID, P_LOCATION_EXTINCTION_PROB, P_RANGE_EXPANSION_PROB, P_HOST_JUMP_PROB, P_SPECIATE_AT_HOST_JUMP
        to_del = []
        for loc in self.geo_range:
            for h in loc.get_hosts_for_para(self):
                if rng.random() < P_LOCATION_EXTINCTION_PROB:
                    to_del.append((loc, h))

        for loc, h in to_del:
            SANITY_CHECK()
            loc.remove_para_from_host_if_present(self, h)
            SANITY_CHECK(False)
        p = find_surrounding(self.geo_range)
        added = []
        for k, v in p.items():
            for src_loc in v:
                src_h_set = src_loc._p2h.get(self, _EMPTY_SET)
                for src_h in src_h_set:
                    if (k in src_h.geo_range) and (rng.random() < P_RANGE_EXPANSION_PROB):
                        added.append((k, src_h))
        if added:
            for loc, src_h in added:
                loc.add_para_for_h_host_range_broken(src_h, [self])
            self.geo_range.update(set([i[0] for i in added]))
        
        new_hosts = {}
        for loc in self.geo_range:
            h_set = set(loc._h2p.keys())
            o_h = loc._p2h[self]
            potential_dest = h_set - o_h
            if potential_dest:
                close_enough_dest_host = [i for i in potential_dest if i.close_enough(o_h)]
                for h in close_enough_dest_host:
                    if rng.random() < P_HOST_JUMP_PROB:
                        new_hosts.setdefault(h, set()).add(loc)
        para_parent = self
        for h, d_geo_range in new_hosts.items():
            if rng.random() < P_SPECIATE_AT_HOST_JUMP:
                para_parent, n = tree._host_jump_speciate(para_parent, para_parent.geo_range, d_geo_range, [h])
            else:
                for loc in d_geo_range:
                    loc.add_para_for_h_host_range_broken(h, [para_parent])
            
        SANITY_CHECK(False)

def _recurse_newick(out, nd, anc, name_pref):
    eob, next_c = nd._get_des_extant_fork()
    if anc:
        br_len = (eob.death_time - anc.death_time)
    if next_c:
        out.write('(')
        for n, c in enumerate(next_c):
            if n != 0:
                out.write(',')
            _recurse_newick(out, c, eob, name_pref)
        if anc:
            out.write('):{b:d}'.format(b=br_len))
        else:
            out.write(')')
    else:
        if anc:
            out.write('{p}{i:d}:{b:d}'.format(p=name_pref, i=eob.index, b=br_len))
        else:
            out.write('{p}{i:d}'.format(p=name_pref, i=eob.index))

class BaseTree(object):
    def newick(self, out, name_pref):
        global CURR_TIME
        for t in self.tips:
            nd = t
            if nd.death_time is None:
                nd.had_none_death_time = True
                nd.death_time = CURR_TIME
            while nd is not None:
                debug('flagging {i:d} as on extant path '.format(i=nd.index))
                nd._on_extant_trace = True
                if nd.index == 0:
                    assert(nd is self.root)
                nd = nd.parent
            assert(self.root._on_extant_trace == True)
        assert(self.root.index == 0)
        assert(self.root._on_extant_trace == True)
        _recurse_newick(out, self.root, None, name_pref)
        out.write(';\n')
        for i in self.tips:
            if getattr(t, 'had_none_death_time', None):
                delattr(t, 'had_none_death_time')
                t.death_time = None

class HostTree(BaseTree):
    def __init__(self, start_time, initial_range, rng):
        self.root = HostLineage(start_time, initial_range)
        self.tips = set([self.root])
        self.max_new_range_dim = 4
        self.new_range_expansion_prob = 0.5 # used only in filling in new ranges. 
        self.prob_dispersal_speciation = 0.1 # Pr(dispersal speciation | speciation)
        self.rng = rng
    def sample_rand_host_loc(self, rng):
        h = rng.sample(self.tips, 1)[0]
        loc = rng.sample(h.geo_range, 1)[0]
        return h, loc

    def go_extinct(self, tip):
        global CURR_TIME
        tip.death_time = CURR_TIME
        self.tips.remove(tip)
        for t in self.tips:
            del t._mrca_times[tip]
        for loc in tip.geo_range:
            try:
                loc.remove_host(tip)
            except: # called in some contexts (during the speciation routine in which the host reassignment has taken place for some locations).
                pass
    def speciate(self, tip):
        global CURR_TIME
        SANITY_CHECK()
        if self.rng.random() < self.prob_dispersal_speciation:
            f, s = self.dispersal_speciate(tip)
        else:
            f, s = self.partition_range_speciate(tip)
        if f is tip:
            SANITY_CHECK()
            return f, s
        self.go_extinct(tip)
        self.tips.add(f)
        self.tips.add(s)
        f._mrca_times[s] = CURR_TIME
        f._mrca_times[f] = CURR_TIME
        SANITY_CHECK()
        return f, s

    def partition_range_speciate(self, tip):
        global CURR_TIME
        nc = len(tip.geo_range)
        if nc < 2:
            return tip, None
        debug('host-partition-range-speciate')
        if Grid.is_tube:
            first_range = set(self.rng.sample(tip.geo_range, nc // 2))
        else:
            if self.rng.random() < 0.5:
                # North/South vicariance
                lat_list = [i.lat for i in tip.geo_range]
                lat_list.sort()
                mid_index = (1 + len(lat_list)) // 2
                mid_lat = lat_list[mid_index]
                below_list, above_list, mid_list = [], [], []
                for loc in tip.geo_range:
                    if loc.lat < mid_lat:
                        below_list.append(loc)
                    elif loc.lat > mid_lat:
                        above_list.append(loc)
                    else:
                        mid_list.append(loc)
            else:
                # East/West vicariance
                lon_list = [i.lon for i in tip.geo_range]
                lon_list.sort()
                mid_index = (1 + len(lon_list)) // 2
                mid_lon = lon_list[mid_index]
                below_list, above_list, mid_list = [], [], []
                for loc in tip.geo_range:
                    if loc.lon < mid_lon:
                        below_list.append(loc)
                    elif loc.lon > mid_lon:
                        above_list.append(loc)
                    else:
                        mid_list.append(loc)
            diff = len(below_list) - len(above_list)
            ms = min(abs(diff), len(mid_list))
            if ms > 0:
                to_transfer = []
                for i in range(ms):
                    to_transfer.append(mid_list.pop())
                if diff < 0:
                    below_list.extend(to_transfer)
                else:
                    above_list.extend(to_transfer)
            to_above, to_below = [], []
            if mid_list:
                to_above = set(self.rng.sample(mid_list, len(mid_list) // 2))
                to_below = list(set(mid_list) - to_above)
                to_above = list(to_above)
            above_list.extend(to_above)
            first_range = set(above_list)
        second_range = tip.geo_range - first_range
        first_daughter = HostLineage(start_time=CURR_TIME,
                                     geo_range=first_range,
                                     parent=tip)
        second_daughter = HostLineage(start_time=CURR_TIME,
                                      geo_range=second_range,
                                      parent=tip)
        for loc in first_range:
            loc.replace_host(tip, first_daughter)
        for loc in second_range:
            loc.replace_host(tip, second_daughter)
        return first_daughter, second_daughter

    def dispersal_speciate(self, tip):
        global CURR_TIME, GRID
        debug('host-dispersal-speciate')
        assert(len(tip.geo_range) > 0)
        same_range_daughter = HostLineage(start_time=CURR_TIME,
                                          geo_range=tip.geo_range,
                                          parent=tip)
        dr = GRID.gen_rand_range(max_ns_extent=self.max_new_range_dim,
                                 max_ew_extent=self.max_new_range_dim,
                                 inc_prob=self.new_range_expansion_prob,
                                 rng=self.rng)
        src_loc = self.rng.sample(tip.geo_range, 1)[0]
        dispersed_parasite_load = src_loc.get_para_for_host(tip)
        for loc in tip.geo_range:
            loc.replace_host(tip, same_range_daughter)
        dispersed_daughter = HostLineage(start_time=CURR_TIME,
                                         geo_range=dr,
                                         parent=tip)
        for loc in dr:
            loc.add_para_for_h_host_range_broken(dispersed_daughter, dispersed_parasite_load)
        return same_range_daughter, dispersed_daughter

class ParasiteTree(BaseTree):
    def __init__(self, start_time, initial_range, host, rng):
        self.root = ParasiteLineage(start_time, initial_range, host)
        self.tips = set([self.root])
        self.rng = rng
    def speciate(self,
                 parasite,
                 h1_range1,
                 h2_range2,
                 geo_speciate_across_hosts=True,
                 other_loc_to_d1=True):
        global GRID, CURR_TIME
        assert(geo_speciate_across_hosts)
        assert(other_loc_to_d1)
        h1, r1 = h1_range1
        h2, r2 = h2_range2
        r = parasite.get_full_range()
        occ_r1 = r - r2
        if (not occ_r1) or not r2:
            return parasite, None
        return self._geo_speciate_across_hosts(parasite, occ_r1, r2)
    def jump_speciate(self, parasite, h_tree):
        h, loc = h_tree.sample_rand_host_loc(self.rng)
        if loc in parasite.geo_range:
            return parasite, None
        parasite.death_time = CURR_TIME
        d1 = ParasiteLineage(CURR_TIME, parasite.geo_range, host=None, parent=parasite)
        d2 = ParasiteLineage(CURR_TIME, set([loc]), host=h, parent=parasite)
        for loc in list(parasite.geo_range):
            h_l = loc.get_hosts_for_para(parasite)
            for h in h_l:
                loc.replace_parasite(parasite, d1)
        self.tips.remove(parasite)
        self.tips.add(d1)
        self.tips.add(d2)
        return d1, d2

    def _geo_speciate_across_hosts(self, parasite, r1, r2):
        parasite.death_time = CURR_TIME
        d1 = ParasiteLineage(CURR_TIME, r1, host=None, parent=parasite)
        d2 = ParasiteLineage(CURR_TIME, r2, host=None, parent=parasite)
        for loc in r1:
            h_l = loc.get_hosts_for_para(parasite)
            for h in h_l:
                loc.replace_parasite(parasite, d1)
        for loc in r2:
            h_l = loc.get_hosts_for_para(parasite)
            for h in h_l:
                loc.replace_parasite(parasite, d2)
        self.tips.remove(parasite)
        self.tips.add(d1)
        self.tips.add(d2)
        return d1, d2

    def _host_jump_speciate(self, parasite, r1, r2, h2):
        parasite.death_time = CURR_TIME
        d1 = ParasiteLineage(CURR_TIME, r1, host=None, parent=parasite)
        d2 = ParasiteLineage(CURR_TIME, r2, host=None, parent=parasite)
        for loc in set(r1):
            h_l = loc.get_hosts_for_para(parasite)
            for h in h_l:
                loc.replace_parasite(parasite, d1)
        for loc in set(r2):
            for h in h2:
                loc.add_para_for_h_host_range_broken(h, [d2])
        self.tips.remove(parasite)
        self.tips.add(d1)
        self.tips.add(d2)
        return d1, d2

    def go_extinct(self, tip):
        global CURR_TIME
        tip.death_time = CURR_TIME
        self.tips.remove(tip)
        for loc in tip.geo_range:
            loc.remove_parasite_broken_range(tip)

    def associations(self, out, hpref, ppref):
        for t in self.tips:
            out.write('{p}{i:d}\t'.format(p=ppref, i=t.index))
            hs = set()
            for loc in t.geo_range:
                hs.update(loc.get_hosts_for_para(t))
            hl = [i.index for i in hs]
            hl.sort()
            for n, h in enumerate(hl):
                if n > 0:
                    out.write(' ')
                out.write('{p}{i:d}'.format(p=hpref, i=h))
            out.write('\n')

def sanity_check(para_list, host_list, non_empty_ranges=True):
    global GRID
    seen_p = set()
    seen_h = set()
    if non_empty_ranges:
        for h in host_list:
            assert(len(h.geo_range) > 0)
        for p in para_list:
            assert(len(p.geo_range) > 0)
    for row in GRID.rows:
        for cell in row:
            hcs = set(cell._h2p.keys())
            pcs = set(cell._p2h.keys())
            pcs_from_h = set()
            hcs_from_p = set()
            for h in hcs:
                seen_h.add(h)
                p_set = cell._h2p[h]
                pcs_from_h.update(p_set)
                assert(cell in h.geo_range)
                for p in p_set:
                    seen_p.add(p)
                    assert(cell in p.geo_range)
                    h_set_from_p = cell._p2h[p]
                    if h not in h_set_from_p:
                        print 'host =', h.__dict__
                        print 'p_set =', p_set
                        print 'p =', p
                        print 'h_set_from_p =', [i.__dict__ for i in h_set_from_p]
                        assert(h in h_set_from_p)
                    hcs_from_p.update(h_set_from_p)
            parasitized_hcs = set([i for i in hcs if cell.get_para_for_host(i)])
            if hcs_from_p != parasitized_hcs:
                print 'cell =', cell.lat, cell.lon
                print 'hcs_from_p =', hcs_from_p
                print 'hcs =', parasitized_hcs
                
                assert(hcs_from_p == parasitized_hcs)
            if pcs_from_h != pcs:
                print 'cell =', cell.lat, cell.lon
                print pcs_from_h
                print pcs
                assert(pcs_from_h == pcs)
            for p in para_list:
                if p not in pcs:
                    assert(cell not in p.geo_range)
            for h in host_list:
                if h not in hcs:
                    assert(cell not in h.geo_range)
    for p in seen_p:
        assert(p in para_list)
    for h in seen_h:
        assert(h in host_list)
    for t1 in host_list:
        for t2 in host_list:
            if t2 is not t1:
                try:
                    eq = (t1._mrca_times[t2] == t2._mrca_times[t1])
                except:
                    eq = False
                if not eq:
                    t1.debug_write(errstream)
                    t2.debug_write(errstream)
                    assert(t1._mrca_times[t2] == t2._mrca_times[t1])

def SANITY_CHECK(non_empty_ranges=True):
    global h_tree, p_tree, SANITY_CHECKS
    if not SANITY_CHECKS:
        return True
    sanity_check(p_tree.tips, h_tree.tips, non_empty_ranges)

def main(h_rng, p_rng, num_hosts):
    global GRID, _HOST_LINEAGE_COUNTER, _PARASITE_LINEAGE_COUNTER, CURR_TIME, h_tree, p_tree, P_RANGE_WIDE_EXTINCTION_PROB
    CURR_TIME = 0
    GRID = Grid(GRID_LENGTH, GRID_LENGTH)
    _HOST_LINEAGE_COUNTER = 0
    _PARASITE_LINEAGE_COUNTER = 0

    init_range = GRID.gen_rand_range(max_ns_extent=max_range_dim,
                                     max_ew_extent=max_range_dim,
                                     inc_prob=.5,
                                     rng=h_rng)
    h_tree = HostTree(start_time=0, initial_range=init_range, rng=h_rng)
    p_tree = ParasiteTree(start_time=0, initial_range=init_range, host=h_tree.root, rng=p_rng)
    if not init_range:
        return h_tree, p_tree
        
    if DEBUGGING_OUTPUT:
        GRID.write_range(errstream, init_range)
    for CURR_TIME in range(1, max_t):
        debug('t={t:d} #h={h:d} #p={p:d}'.format(
                    t=CURR_TIME,
                    h=len(h_tree.tips),
                    p=len(p_tree.tips)))
        to_speciate = set()
        print_ranges = DEBUGGING_OUTPUT and (len(h_tree.tips) == 1)
        if print_ranges:
            for t in h_tree.tips:
                errstream.write('Before speciation\n')
                GRID.write_range(errstream, t.geo_range)
        ####
        # Host speciation...
        for h_l in h_tree.tips:
            if h_rng.random() < H_SPECIATION_PROB:
                to_speciate.add(h_l)
        # It is not great to stop at the first time
        #   you were going to exceed the number of hosts, 
        #   but this should work for our purposes...
        if (len(to_speciate) > 0) and (len(h_tree.tips) + len(to_speciate) > num_hosts):
            if len(h_tree.tips) == num_hosts:
                break
            if len(to_speciate) > 1:
                n_events = num_hosts - len(h_tree.tips)
                to_speciate = set(list(to_speciate)[:n_events])
            
        for h_l in to_speciate:
            assert(len(h_l.geo_range) > 0)
            d1, d2 = h_tree.speciate(h_l)
            if print_ranges:
                for t in h_tree.tips:
                    errstream.write('After speciation\n')
                    GRID.write_range(errstream, t.geo_range)
                if (raw_input('Continue? ').lower() == 'n'):
                    sys.exit(1)
            SANITY_CHECK()
            if d2 is not None: # speciation happened (not just tried)
                if H_PROB_PARA_SP_GIVEN_HOST_SP > 0.0:
                    d1_para = d1.collect_para2range()
                    d2_para = d2.collect_para2range()
                    ak = set(d1_para.keys())
                    ak.update(set(d2_para.keys()))
                    for para_lineage in ak:
                        p1 = d1_para.get(para_lineage)
                        p2 = d2_para.get(para_lineage)
                        if (p1 is not None) and (p2 is not None):
                            p_tree.speciate(para_lineage, (d1, p1), (d2, p2))
                            SANITY_CHECK()
        
        #####   
        # Host range changes...
        gone_extinct = []
        
        for t in h_tree.tips:
            t.disperse(h_tree.rng)
            if (not t.geo_range) and (t.death_time is None):
                gone_extinct.append(t)
        for t in gone_extinct:
            h_tree.go_extinct(t)
        if not h_tree.tips:
            break
        SANITY_CHECK()
        
        gone_extinct = []
        old_p_tips = list(p_tree.tips) # make a copy, because disperse can lead to speciation...
        for p in old_p_tips:
            if p_tree.rng.random() < P_RANGE_WIDE_EXTINCTION_PROB:
                gone_extinct.append(p)
            else:
                p.disperse(p_tree.rng, p_tree)
                if (not p.geo_range) and (p.death_time is None):
                    gone_extinct.append(p)
                if P_MAX_RANGE_IN_EFFECT:
                    if len(p.geo_range) > P_MAX_RANGE:
                        num_to_del = len(p.geo_range) - P_MAX_RANGE
                        to_del =  p_tree.rng.sample(p.geo_range, num_to_del)
                        for loc in to_del:
                            loc.remove_parasite_broken_range(p)
                            p.geo_range.remove(loc)
        for p in gone_extinct:
            p_tree.go_extinct(p)
        if not p_tree.tips:
            break
        if P_LONG_DIST_PROB > 0.0:
            tl = set(p_tree.tips)
            for p in tl:
                if p_tree.rng.random() < P_LONG_DIST_PROB:
                    if p_tree.rng.random() < P_SP_ON_LONG_DIST_PROB:
                        p_tree.jump_speciate(p, h_tree)
                    else:
                        h, loc = h_tree.sample_rand_host_loc(p_tree.rng)
                        loc.add_para_for_h_host_range_broken(h, [p])
        if len(p_tree.tips) > (ABORT_P_TIPS_MULTIPLIER * MAX_NUM_PARA):
            break
        SANITY_CHECK()
    return h_tree, p_tree

if __name__ == '__main__':
    from random import Random
    import sys
    import codecs
    import os
    try:
        from ConfigParser import SafeConfigParser
    except:
        from configparser import SafeConfigParser
    
    errstream = codecs.getwriter('utf-8')(sys.stderr)
    try:
        h_seed = int(sys.argv[1])
        p_seed = int(sys.argv[2])    
        NUM_HOSTS = int(sys.argv[3])
        MIN_NUM_PARA = int(sys.argv[4])
        MAX_NUM_PARA = int(sys.argv[5])
        # if the number of parasite species exceeds ABORT_P_TIPS_MULTIPLIER*MAX_NUM_PARA
        #   that simulation rep will be aborted to save on computation on a replicate
        #   that is very unlikely to produce a replicate with the constraints on the
        #   number of parasite species.
        ABORT_P_TIPS_MULTIPLIER = 10 
        if MIN_NUM_PARA > MAX_NUM_PARA:
            error('minimum # of parasites cannot be larger than maximum # parasites')
        assert(MIN_NUM_PARA <= MAX_NUM_PARA)
        config_fp = sys.argv[6]
        if not os.path.exists(config_fp):
            error('Config file {f} does not exist'.format(f=config_fp))
        assert(os.path.exists(config_fp))
        parser = SafeConfigParser()
        parser.read(config_fp)
        # Read global section
        # Tube is GRID_LENGTH on each side, and GRID_LENGTH top to bottom
        GRID_LENGTH = parser.getint('global', 'grid-length')
        assert(GRID_LENGTH > 0)
        max_t = parser.getint('global', 'max-number-iterations')
        assert(max_t > 0)
        # Read host section
        max_range_dim = parser.getint('host', 'max-init-range-width')
        assert(max_range_dim > 0)
        assert(max_range_dim <= GRID_LENGTH)
        H_SPECIATION_PROB = parser.getfloat('host', 'speciation-prob')
        assert(H_SPECIATION_PROB > 0.0)
        H_LOCATION_EXTINCTION_PROB = parser.getfloat('host', 'local-extinction-prob')
        assert(H_LOCATION_EXTINCTION_PROB >= 0.0)
        H_RANGE_EXPANSION_PROB = parser.getfloat('host', 'range-expansion-prob')
        assert(H_RANGE_EXPANSION_PROB >= 0.0)
        # Read parasite section
        H_PROB_PARA_SP_GIVEN_HOST_SP = parser.getfloat('parasite', 'speciation-prob-given-host-speciation')
        assert(H_PROB_PARA_SP_GIVEN_HOST_SP >= 0.0)
        P_LOCATION_EXTINCTION_PROB = parser.getfloat('parasite', 'local-extinction-prob')
        assert(P_LOCATION_EXTINCTION_PROB >= 0.0)
        P_RANGE_EXPANSION_PROB = parser.getfloat('parasite', 'range-expansion-prob')
        assert(P_RANGE_EXPANSION_PROB >= 0.0)
        P_HOST_JUMP_PROB = parser.getfloat('parasite', 'new-host-infection-prob')
        assert(P_HOST_JUMP_PROB >= 0.0)
        PHYLOGENETIC_HOST_JUMPS = parser.getboolean('parasite', 'new-host-depends-on-phylogeny')
        host_infection_patristic_threshold = parser.getfloat('parasite', 'host-infection-patristic-threshold')
        P_SPECIATE_AT_HOST_JUMP = parser.getfloat('parasite', 'speciate-on-new-host-infection-prob')
        assert(P_SPECIATE_AT_HOST_JUMP >= 0.0)
        P_RANGE_WIDE_EXTINCTION_PROB = parser.getfloat('parasite', 'range-wide-extinction-prob')
        assert(P_RANGE_WIDE_EXTINCTION_PROB >= 0.0)
        P_MAX_RANGE = parser.getint('parasite', 'max-geo-range')
        assert(P_MAX_RANGE > 0)
        P_MAX_RANGE_IN_EFFECT = (P_MAX_RANGE < (GRID_LENGTH*GRID_LENGTH))
        P_LONG_DIST_PROB = parser.getfloat('parasite', 'long-dist-colonize-prob')
        assert(P_LONG_DIST_PROB >= 0.0)
        P_SP_ON_LONG_DIST_PROB = parser.getfloat('parasite', 'speciate-on-long-dist-colonize-prob')
        assert(P_SP_ON_LONG_DIST_PROB >= 0.0)

        

    except Exception, x:
        import traceback
        xs = traceback.format_exc()
        cfg_help='''Config file:
    The configfile should have the standard python config file syntax. 
    See the example at the top of
        http://docs.python.org/release/2.7/library/configparser.html 
    for details. For this program, the config should have:
    A "global" section with the following fields:
        "grid-length" int The width and height of the grid.
    A "host" section with the following fields:
        "speciation-prob" float Per time step probability of speciation for
            each host lineage.
            
    
    See example.cfg for and example.
'''
        sys.exit('''
    
Program:
    sim-host-parasite.py  Copyright (C) 2014 Mark T. Holder mtholder@gmail.com
    This program comes with ABSOLUTELY NO WARRANTY. See LICENSE.txt for details.

    Simple simulator host and parasite genealogies on a geographic grid

Prerequisites:
    None other than python 2 standard library. Tested on python 2.7

Usage:

    python host-seed parasite-seed #-hosts min-#-para-sp max-#-para-sp configfile

  to simulate a host tree of #-hosts species using a rng with host-seed. A
  parasite tree will be simulated using parasite-seed. The simulations will be
  repeated until a realization produces a parasite tree with a number of species
  between min-#-para-sp and max-#-para-sp using the settings in configfile
  to parameterize the details of the simulation
  
  For example:
    python sim-host-parasite.py $RANDOM $RANDOM 50 40 45 example.cfg
  will produces a host tree of 50 species and a parasite tree with a number of
  tips in the interval [40, 45]
  
  Warning messages are printed to stderr.
  The output is printed to stdout.

''' + cfg_help + '\n' + 'Error. Exception: ' + xs + '\n')
    print h_seed, p_seed
    DEBUGGING_OUTPUT = os.environ.get('DEBUGGING_OUTPUT', '0') != '0'
    SANITY_CHECKS = os.environ.get('SANITY_CHECKS', '0') != '0'
    h_rng = Random(h_seed)
    p_rng = Random(p_seed)
    while True:
        h_tree, p_tree = main(h_rng,
                              p_rng,
                              num_hosts=NUM_HOSTS)
        if len(h_tree.tips) == NUM_HOSTS:
            num_para = len(p_tree.tips)
            if num_para < MIN_NUM_PARA:
                warn('#para = {p:d} < {t:d}'.format(p=num_para, t=MIN_NUM_PARA))
            elif num_para > MAX_NUM_PARA:
                warn('#para = {p:d} > {t:d}'.format(p=num_para, t=MAX_NUM_PARA))
            else:
                break
    out = sys.stdout
    h_tree.newick(out, 'h')
    p_tree.newick(out, 'p')
    p_tree.associations(out, 'h', 'p')

# Simple simulator host and parasite genealogies on a geographic grid
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
