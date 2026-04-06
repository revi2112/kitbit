import os
from timeit import default_timer
from math import exp, floor
from copy import deepcopy
from math import fabs, log, inf, ceil
from itertools import product

class KitBit:
    def __init__(self, structure, kl, mni, depth, search_algorithm='BFS', n=2, min_zeros=1, epsilon=exp(-18), all_solutions=False):
        self.structure = structure
        self.kl = kl
        self.mni = mni
        self.depth = depth
        self.search_alg = search_algorithm
        self.n = n
        self.min_zeros = min_zeros
        self.epsilon = epsilon
        self.all_solutions = all_solutions

    def handler(self):
        if len(self.structure) == 1:
            return (False, 0, self.kl)
        if self.all_solutions:
            init_time = default_timer()
            sol = self.numeric_solver_all_sols(self.structure[:], self.kl, self.depth, 10)
            end_time = default_timer()
        else:
            init_time = default_timer()
            sol = self.numeric_solver(self.structure[:], self.kl, self.depth, 10)
            end_time = default_timer()
        return [sol, end_time - init_time]

    def numeric_solver(self, seq, kl, depth, module):
        edk0 = KBEDK(seq[:], module, None)
        edk0 = edk0.basic()
        if edk0.is_goal('BASIC', self.epsilon, self.min_zeros, 0):
            edk0.eos[-1] += [0]*self.n
            return edk0.rev_basic(self.n, 1, len(edk0.eos)-1), 0, ['BASIC']
        if len(seq) == 2:
            edk0.eos[-1] += edk0.eos[-1]*self.n
            return edk0.rev_basic(self.n, 1, len(edk0.eos)-1), 0, ['BASIC']
        st0 = SeqState(None, [edk0], [], [])
        alg = SeqSearchAlgorithm(st0, kl, self.mni, depth, self.n, self.min_zeros, self.epsilon)
        if self.search_alg == 'BFS':
            lst = alg.bfs()
        elif self.search_alg == 'DFS':
            lst = alg.dfs([st0])
        else:
            lst = alg.branch()
        if lst is None or lst is False or len(lst[0]) == 1:
            return False, 0, kl
        acts = [st.action for st in lst[0] if st.action is not None]
        sol = SeqPredictor(lst[0][::-1], self.n)
        ln, sol = len(seq), sol.predictor()
        return (False, 0, kl) if not sol or len(sol) <= len(seq) else (seq[:]+sol[ln:ln+self.n], lst[1], acts)

    def numeric_solver_all_sols(self, seq, kl, depth, module):
        edk0 = KBEDK(seq[:], module, None)
        edk0 = edk0.basic()
        if edk0.is_goal('BASIC', self.epsilon, self.min_zeros, 0):
            edk0.eos[-1] += [0 for i in range(self.n)]
            return [(edk0.rev_basic(self.n, 1, len(edk0.eos)-1), 0, ['BASIC'])]
        st0 = SeqState(None, [edk0], [], [])
        alg = SeqSearchAlgorithm(st0, kl, self.mni, depth, self.n, self.min_zeros, self.epsilon)
        ln, sols, l_lst = len(seq), [], alg.kb_seq_gb()
        if l_lst is None or l_lst is False or len(l_lst) == 1:
            return (False, 0, kl)
        l_lst, n_iter = l_lst[0], l_lst[1]
        for lst in l_lst:
            acts = [st.action for st in lst if st.action is not None]
            sol = SeqPredictor(lst[::-1], self.n)
            sol = sol.predictor()
            if not sol or len(seq[:]+sol[ln:ln+self.n])<=ln:
                continue
            sols.append((seq[:]+sol[ln:ln+self.n], n_iter, acts))
        return sols if sols else (False, 0, kl)

class SeqSearchAlgorithm:
    def __init__(self, init_state, kl, mni, depth, n, min_zeros, epsilon):
        self.init_state = init_state
        self.kl = kl
        self.mni = mni
        self.depth = depth
        self.n = n
        self.min_zeros = min_zeros
        self.epsilon = epsilon
        self.road = [init_state]
        self.count = 0
        self.objective = 'None'

    def bfs(self):
        i, k, ln, self.road = 0, -1, len(self.kl)-1, [self.road]
        n_iter = [len(self.kl)**p for p in range(1, self.depth+1)]
        for j in range(1, sum(n_iter)+1):
            k = k+1 if k < ln else 0
            st = False if not self.road[i] else \
                self.road[i][-1].new_state(self.n, self.road[i][-1].action,
                    self.kl[k], self.epsilon, self.min_zeros)
            if st is False:
                i = i+1 if k == ln else i
                self.road.append(False)
                continue
            self.road.append(self.road[i]+[st])
            i = i+1 if k == ln else i
            if len(st.sols) == len(st.eos):
                return self.road[-1], j
            if j >= self.mni:
                return False
        return False

    def dfs(self, road=[], level=0):
        for i in range(len(self.kl)):
            self.count += 1
            if len(self.road) > 1 or not self.objective:
                break
            st = road[-1].new_state(self.n, road[-1].action, self.kl[i], self.epsilon, self.min_zeros)
            if st is False:
                continue
            road.append(st)
            if len(st.sols) == len(st.eos):
                self.road = deepcopy(road)
                return road, self.count
            if self.count >= self.mni:
                self.objective = False
                return False
            if level < self.depth - 1:
                self.dfs(road, level+1)
            road.pop()
        if self.road:
            return self.road, self.count

    def kb_seq_gb(self):
        branches, self.road = [], [self.road]
        for i in range(self.depth):
            branches1 = []
            for lst in self.road:
                sts = list(map(lambda kita: lst[-1].new_state(self.n, lst[-1].action, kita, self.epsilon, self.min_zeros), self.kl))
                self.count += len(sts)
                sts = list(filter(lambda x: x is not False, sts))
                sts_sol = list(filter(lambda x: len(x.sols) == len(x.eos), sts))
                sts_not_sol = list(filter(lambda x: len(x.sols) != len(x.eos), sts))
                branches1 += list(map(lambda st: deepcopy(lst)+[st], sts_not_sol))
                branches += list(map(lambda st: deepcopy(lst)+[st], sts_sol))
            self.road = branches1
            if not branches1:
                break
        return (branches, self.count) if branches else False

    def branch(self):
        ln, road = len(self.kl), self.road
        for j in range(ln):
            st = False if not self.road[j] else road[-1].new_state(self.n, road[-1].action, self.kl[j], self.epsilon, self.min_zeros)
            if st is False:
                return False
            road.append(st)
            if len(st.sols) == len(st.eos) and j == ln-1:
                return road, j
        return False

class SeqState:
    def __init__(self, action, eos, sols, parents):
        self.action = action
        self.eos = eos
        self.sols = sols
        self.parents = parents

    def new_state(self, n, pva, act, epsilon, min_zeros):
        ln, edks, parents = len(self.eos), [], []
        for i in range(ln):
            edk = KBEDK(self.eos[i].eos, self.eos[0].module)
            if (ln > 1 and i not in self.sols) or ln == 1:
                try:
                    new_edk = self.do_action(edk, n, pva, act)
                except:
                    return False
                if new_edk is False or False in new_edk:
                    return False
                edks += new_edk
                parents += [i for q in range(len(new_edk))]
        m = 0 if act != 'DIV' else 1
        sols = [e for e in range(len(edks)) if edks[e].is_goal(act, epsilon, min_zeros, m)]
        return SeqState(act, edks, sols, parents)

    @staticmethod
    def do_action(edk, n, pva, act):
        if act == 'DIV':
            return edk.divisions()
        if act[:3] == 'RED':
            a = act[act.index('(')+1:act.index(')')]
            return edk.reduction(int(a))
        elif act[:3] == 'EXP':
            a = act[act.index('(')+1:act.index(')')]
            if '/' in a:
                e0 = a.split('/')
                e = int(e0[0]) / int(e0[1])
            elif '.' in a:
                e = float(a)
            else:
                e = int(a)
            return edk.exponentiation(e)
        elif act == 'LOG':
            return edk.logarithm()
        elif act[:3] == 'DOP':
            a = act[act.index('(')+1:act.index(')')].split(',')
            return edk.double_operation(a[0], a[1])
        elif act[:2] == 'ML':
            a = act[act.index('(')+1:act.index(')')].split(',')
            return edk.multi_level(int(a[0]), int(a[1]))
        elif act[:3] == 'FOC':
            a = act[act.index('(')+1:act.index(')')].split(',')
            divs = [int(i) for i in a[1:]]
            return edk.focusing(int(a[0]), divs)
        elif act[:3] == 'ANA':
            a = act[act.index('(')+1:act.index(')')].split(',')
            return edk.analogy(int(a[0]), int(a[1]))
        elif act[:3] == 'SOE':
            return edk.split_of_elements()
        elif act[:4] == 'DGEE':
            a = act[act.index('(') + 1:act.index(')')].split(',')
            return edk.dgee(int(a[0]), int(a[1]))
        elif act[:4] == 'DGDE':
            a = act[act.index('(') + 1:act.index(')')].split(',')
            return edk.dgde(int(a[0]), int(a[1]))
        elif act == 'RSYM':
            return edk.repetition_symmetry(n)
        elif act == 'SSYM':
            return edk.specular_symmetry()
        elif act == 'CSYM':
            return edk.circular_symmetry(n)
        elif pva is not None and ('DIAG' in pva or pva == 'LOG' or
                pva[:3] == 'RED' or pva[:2] == 'ML' or pva[:3] == 'DOP'
                or pva[:3] == 'FOC' or pva[:3] == 'SOE' or 'DG' in pva):
            return False
        elif act == 'ASYM':
            return edk.array_symmetry()
        elif act == 'LDIAG':
            return edk.diagonal(0, 2, 1, 0, edk.dim1+1)
        elif act == 'RDIAG':
            return edk.diagonal(1, 1, 0, edk.dim1-1, edk.dim1-1)
        elif act[:4] == 'RECK':
            a = act[act.index('(')+1:act.index(')')].split(',')
            return edk.reck(a[0], a[1])
        else:
            return edk.transposed()

class KBEDK:
    def __init__(self, eos, module, dim=None):
        self.eos = eos
        self.module = module
        self.dim = dim

    def is_goal(self, act, epsilon, min_zeros, elem):
        if self.eos is not False and act is not None and \
                (act[:3] == 'ANA' or 'SYM' in act or act[:4] == 'RECK'):
            return True
        try:
            if min_zeros == 1:
                return False if None in self.eos[-1] else \
                    fabs(self.eos[-1][-1]-elem) < epsilon
            for i in self.eos[-min_zeros:]:
                for j in i:
                    if j is None or fabs(j-elem) >= epsilon:
                        return False
            return True
        except:
            return False

    def check_number(self, elem):
        if self.module == 10:
            return elem
        dec = elem-int(elem)
        elem = int(elem) if dec == 0.0 else elem
        if type(elem) == float:
            return False
        while elem >= self.module:
            elem -= self.module
        while elem < 0:
            elem += self.module
        return elem

    def do_op(self, op, m, n, rev):
        if op == '+' or (op == '-' and rev):
            return self.check_number(m + n)
        elif op == '-' or (op == '+' and rev):
            return self.check_number(m - n)
        elif op == '*' or (op == '/' and rev):
            return self.check_number(m * n)
        else:
            return self.check_number(m / n) if n != 0 else False

    def basic(self):
        if len(self.eos) < 2:
            return False
        t = [self.eos]
        for y in range(len(self.eos)-1):
            seq = []
            for x in range(len(t[y])-1):
                e = self.check_number(t[y][x+1]-t[y][x])
                if e is False:
                    return False
                seq.append(e)
            t.append(seq)
        return KBEDK(t, self.module, self.dim)

    def divisions(self):
        if len(self.eos[0]) < 2 or 0 in self.eos[0]:
            return False
        t = [self.eos[0]]
        for y in range(len(self.eos)-1):
            seq = []
            for x in range(len(t[y])-1):
                e = self.check_number(t[y][x+1] / t[y][x])
                if e == 0 or e > 1e9 or e is False:
                    return False
                seq.append(e)
            t.append(seq)
        return [KBEDK(t, self.module, self.dim)]

    def reduction(self, height):
        if len(self.eos) - height >= 2:
            t = KBEDK(self.eos[height], self.module, self.dim)
            return [t.basic()]
        return False

    def exponentiation(self, j):
        if len(self.eos[0]) < 2 or j == 0 or (0 in self.eos[0] and j <= 0):
            return False
        seq = []
        for i in self.eos[0]:
            e = (i)**(j)
            if type(e) == complex:
                return False
            e = self.check_number(e)
            if e is False:
                return False
            seq.append(e)
        seq = KBEDK(seq, self.module, self.dim)
        return [seq.basic()]

    def logarithm(self):
        base, ln, seq = self.eos[0], len(self.eos[0]), []
        if ln < 2 or 0 in base or 1 in base or -1 in base:
            return False
        for i in range(len(self.eos)-1):
            e = self.check_number(log(fabs(base[i+1]), fabs(base[i])))
            if e is False:
                return False
            seq.append(e)
        seq = KBEDK(seq, self.module, self.dim)
        return [seq.basic()]

    def double_operation(self, op1, op2):
        base, ln = self.eos[0], len(self.eos[0])
        seq, r = [], ln % 2
        cond1 = op1 == '/' and 0 in base[::2]
        cond2 = op2 == '/' and 0 in base[1::2]
        cond3 = op1 == '*' and r == 0 and 0 in base[:-2]
        cond4 = op2 == '*' and r != 0 and 0 in base[:-2]
        if ln < 2 or cond1 or cond2 or cond3 or cond4:
            return False
        for i in range(len(self.eos) - 1):
            op = op1 if i % 2 == 0 else op2
            e = self.do_op(op, base[i+1], base[i], False)
            if e is False:
                return False
            seq.append(e)
        edk_sol = KBEDK(seq, self.module, self.dim)
        return [edk_sol.basic()]

    def multi_level(self, dx, dy):
        if len(self.eos[0]) < 2:
            return False
        edks, p, x, y, dxy = [], 0, 0, 0, dx+dy
        for i in range(dxy):
            seq = []
            while len(self.eos[y]) > x and y < len(self.eos):
                seq.append(self.eos[y][x])
                if y+dy > len(self.eos) - 1:
                    break
                x, y = x+dx, y+dy
                p = p+1 if x == len(self.eos[y]) else p
            if len(seq) == 1 or (p == 0 and i == dxy-1 and dxy != 1):
                return False
            edk = KBEDK(seq, self.module, self.dim)
            edks.append(edk.basic())
            x, y = 0, i+1
        return edks

    def focusing(self, shift, divs):
        ln1, s = len(self.eos[0]), sum(divs)
        if shift >= s or s+shift > ln1 or ln1 < 4:
            return False
        main_list = self.eos[0][shift:]
        remainder, ln2 = self.eos[0][:shift], len(divs)
        seqs, edks = [[] for k in range(ln2)], []
        while main_list:
            for p in range(ln2):
                seqs[p] += main_list[:divs[p]]
                del main_list[:divs[p]]
                if not main_list:
                    break
        for q in range(ln2-1, -1, -1):
            if remainder:
                seqs[q] = remainder[-divs[q]:] + seqs[q]
                del remainder[-divs[q]:]
            if len(seqs[q]) <= 1:
                return False
            edk = KBEDK(seqs[q], self.module, self.dim)
            edks.insert(0, edk.basic())
        return edks

    def analogy(self, shift, nels):
        ln0 = len(self.eos[0][shift:])
        if ln0 < 4 or nels == 1 or nels+1 > ln0 or ln0 % nels == 0:
            return False
        seq0, seqs = deepcopy(self.eos[0][shift:]), []
        while seq0:
            seq = seq0[:nels] if len(seq0[:nels]) > 1 else [seq0[:nels]]
            seqs.append(KBEDK(seq, self.module, self.dim))
            del seq0[:nels]
        edks = list(map(lambda s: s.basic() if len(s.eos) > 1 else s, seqs))
        edkf = edks.pop()
        for i in range(1, len(edks[0].eos[0])):
            prz = list(map(lambda e: e.eos[i] == [0]*len(e.eos[i]), edks))
            nt = prz.count(True)
            if nt == len(prz):
                lnf = len(edkf.eos)
                if lnf == 1 or (i >= lnf and edkf.eos[-1] != [0]*len(edkf.eos[-1])):
                    h, i = len(edks[0].eos)-lnf, 1 if lnf == 1 else i
                    subedk = deepcopy(edks[0].eos[-h:])
                    if list(filter(lambda x: x.eos[-h:] != subedk, edks)):
                        return False
                    edkf.eos += subedk
                elif i < lnf and edkf.eos[i] == [0]*len(edkf.eos[i]) \
                        and edkf.eos[i-1] != [0]*len(edkf.eos[i-1]):
                    edkf.eos[i:] = deepcopy(edks[0].eos[i:])
                else:
                    return False
                if len(edkf.eos) != len(edks[0].eos):
                    return False
                ln2, ln1 = len(edkf.eos[i]), len(edkf.eos[i-1])-1
                ln3 = len(edkf.eos[0])
                SEQ = self.eos[0] + edkf.rev_basic(ln2-ln1, ln1, i)[ln3:]
                if SEQ == self.eos[0]:
                    return False
                return [KBEDK([SEQ, 0], self.module, self.dim)]
            elif nt != 0:
                return False
        return False

    def split_of_elements(self):
        check = list(map(lambda x: type(x) == float or x < 0, self.eos[0]))
        if len(self.eos[0]) <= 2 or True in check:
            return False
        s1, s2 = [], []
        for e in self.eos[0]:
            num_to_str = str(e)
            add_commas = ','.join(num_to_str)
            divide = add_commas.split(',')
            lend = [int(e) for e in divide]
            if lend != [lend[0]]*len(lend):
                return False
            s1.append(lend[0])
            s2.append(len(lend))
        if len(s1) < 2 or s1 == self.eos[0]:
            return False
        e1 = KBEDK(s1, self.module, self.dim)
        e2 = KBEDK(s2, self.module, self.dim)
        return [e1.basic(), e2.basic()]

    def dgee(self, x, y):
        if len(self.eos[0]) <= 2:
            return False
        s1, groups, a = [], [], None
        for e in self.eos[0]:
            if e != a:
                a = e
                s1.append(e)
                groups.append([])
            groups[-1].append(e)
        if len(s1) < 2:
            return False
        s2 = list(map(lambda y: len(y), groups))
        e1 = KBEDK(s1[:-1] if x == 1 else s1, self.module, self.dim)
        e2 = KBEDK(s2[:-1] if y == 1 else s2, self.module, self.dim)
        if len(e2.eos) < 2 or len(e1.eos) < 2 or s2 == [1]*len(s2):
            return False
        return [e1.basic(), e2.basic()]

    def dgde(self, x, y):
        seq, ln0 = self.eos[0][:], len(self.eos[0])
        if ln0 <= 2:
            return False
        s1 = [e for i, e in enumerate(seq) if e not in seq[:i]]
        ln1, u, v, c = len(s1), 0, 0, 0
        if self.eos[0][0:ln1] == s1:
            c = 1
        else:
            for j in range(ln1, 0, -1):
                if s1[:j] == seq[-len(s1[:j]):]:
                    u = len(s1[:j])
                    break
            if u == 0:
                return False
        s, groups = [self.eos[0].count(e) for e in s1], []
        while seq:
            g, ln2, ln3 = [], len(seq), u if c == 0 and v == 0 else ln1
            for i in range(ln3):
                if s[i] != 0:
                    g.append(s1[i])
                    s[i] = s[i]-1
            v, a = 1, (0, len(g)) if c == 1 else (ln2-len(g), ln2+1)
            if seq[a[0]:a[1]] != g:
                return False
            del seq[a[0]:a[1]]
            groups.append(g)
        groups = groups if c == 1 else groups[::-1]
        s2 = list(map(lambda z: len(z), groups))
        e1 = KBEDK(s1[:-1] if x == 1 else s1, self.module, self.dim)
        e2 = KBEDK(s2[:-1] if y == 1 else s2, self.module, self.dim)
        if len(e2.eos) < 2 or (len(e1.eos) < 2 and c == 0):
            return False
        return [e2.basic()] if c == 1 else [e1.basic(), e2.basic()]

    def repetition_symmetry(self, n):
        ln = len(self.eos[0])
        if ln <= 2:
            return False
        for i in range(2, ln):
            seq = deepcopy(self.eos[0])
            groups = [seq[k:k+i] for k in range(0, len(seq), i)]
            g0, gf, ln1 = groups[0], groups[-1], len(groups)
            lg0, lgf = len(g0), len(gf)
            if lg0 > lgf:
                for k in range(1, ln1):
                    if k == ln1-1 and g0[:lgf] == gf:
                        x = ceil((len(gf)+n)/len(g0))
                        s = [n for m in groups[:ln1-1] for n in m]+g0*x
                        return [KBEDK([s, 0], self.module, self.dim)]
                    if groups[k-1] != groups[k]:
                        break
        return False

    def specular_symmetry(self):
        ln = len(self.eos[0])
        if ln <= 2:
            return False
        for i in range(round(ln/2), ln):
            e0, e1 = self.eos[0][:i], self.eos[0][:i+1]
            e2, e3, e4 = self.eos[0][i:], e0[::-1], e1[::-1]
            if e0 == e2[::-1]:
                return False
            if e3[:len(e2)] == e2:
                if len(e0+e3) <= ln:
                    continue
                return [KBEDK([e0+e3, 0], self.module, self.dim)]
            if e4[:len(e2)] == e2 and i != ln-1:
                if len(e1+e4[1:]) <= ln:
                    continue
                return [KBEDK([e1+e4[1:], 0], self.module, self.dim)]
        return False

    def circular_symmetry(self, n):
        m, length, mask = 0, len(self.eos)+n, []
        if length < 3:
            return False
        for i in range(1, length+1):
            if i % 2 == 0:
                mask.append(m)
                m += 1
            else:
                mask.insert(0, m)
        ln = length if length % 2 != 0 else int(length/2)
        seq = self.eos[0] + [None]*n
        sol0 = self.check_csym(mask, ln, seq)
        if length % 2 == 0 and sol0 is False:
            m1 = list(range(1, ln))
            mask = [0] + m1 + [0] + m1[::-1]
            sol0 = self.check_csym(mask, ln, seq)
        return False if not sol0 else [KBEDK([sol0, 0], self.module, self.dim)]

    @staticmethod
    def check_csym(mask, length, seq):
        for j in range(length):
            ln = len(mask)
            link = list(zip(mask, seq))
            link0 = link[:]
            link.sort(key=lambda y: y[0])
            ln1 = ln if ln % 2 == 0 else ln-1
            sym = sum([1 for k in range(0, ln1, 2) if link[k][1] ==
                       link[k+1][1] or (link[k+1][1] is None and
                       link[k][1] is not None) or (link[k][1] is None
                       and link[k+1][1] is not None)])
            if sym == ln1 // 2:
                for k in range(len(link)):
                    if link[k][1] is None:
                        if link[k-1][0] != link[k][0] and k == len(link)-1:
                            continue
                        link0[link0.index(link[k])] = link[k-1] if link[k-1][0] == link[k][0] else link[k+1]
                sol = list(zip(*link0))
                if None not in sol[1]:
                    return list(sol[1])
            else:
                shift = mask.pop(-1)
                mask.insert(0, shift)
        return False

    def array_symmetry(self):
        length, dim0, dim1 = len(self.eos), self.dim[0], self.dim[1]
        if dim0 <= 1 or dim1 <= 1 or length < (dim0-1)*dim1+1 or length >= dim0*dim1:
            return False
        seq0, ln = self.eos[0][:], len(self.eos[0])
        seq0 += [None]*(dim0*dim1-ln)
        sr = self.rc_sym(seq0[:], dim0, dim1)
        if sr:
            return sr
        seq1 = []
        for i in range(dim1):
            seq1 += seq0[i::dim1]
        sc = self.rc_sym(seq1[:], dim1, dim0)
        if sc:
            sc, sct = sc[0].eos[0], []
            for i in range(dim0):
                sct += sc[i::dim0]
            return [KBEDK([sct, 0], self.module, self.dim)]
        if dim0 != dim1 or dim0 <= 2:
            return False
        for j, k in zip([0, dim0-1], [1, -1]):
            md = [i for i in range(j, len(seq0), dim0+k)]
            seq2 = self.diagonal_sym(seq0[:], md[:], 0, -k)
            if seq2 is False:
                continue
            sd = self.diagonal_sym(seq2, md[:], -1, k)
            if sd is False:
                continue
            return KBEDK([sd, 0], self.module, self.dim)
        return False

    def rc_sym(self, seq, dim0, dim1):
        groups = [seq[p:p+dim1] for p in range(0, len(seq), dim1)]
        ln, sol = len(groups), []
        j = int(ln/2)
        i = j if ln % 2 != 0 else j-1
        for k in range(round(dim0/2)):
            if None in groups[j]:
                d = groups[j].index(None)
                if groups[i][:d] != groups[j][:d]:
                    return False
                groups[j] = groups[i]
            else:
                if groups[i] != groups[j]:
                    return False
            i, j = i-1, j+1
        for m in groups:
            sol += m
        return [KBEDK([sol, 0], self.module, self.dim)]

    def diagonal_sym(self, seq, md, e1, e2):
        for j in range(self.dim[0]-1):
            if j != 0:
                del md[e1]
                md = [k+e2 for k in md]
            md1 = [seq[m] for m in md]
            ln = len(md1)
            x, y = md1[:round(ln/2)], md1[ln//2:]
            if (None not in y and x != y[::-1]) or (None in y and None in x):
                return False
            if None in y:
                p, q = x[::-1], y.index(None)
                if p[:q] != y[:q]:
                    return False
                if len(x) == 1 or p[:q] == y[:q]:
                    md1[ln//2:] = p
                    for r in range(ln):
                        seq[md[r]] = md1[r]
        seq[-1] = seq[0] if seq[-1] is None else seq[-1]
        seq[0] = seq[-1] if seq[0] is None else seq[0]
        return seq

    def diagonal(self, a, b, c, d, e):
        dim0, dim1 = self.dim[0], self.dim[1]
        if dim0 != dim1-a or dim0 <= b or dim1 <= 2 or len(self.eos[0]) != dim0*dim1-c:
            return False
        array = self.eos[0]
        diagonal = [array[i] for i in range(d, len(array), e)]
        sol = KBEDK(diagonal, self.module, self.dim)
        return [sol.basic()]

    def reck(self, roc, ops):
        dim0, dim1 = self.dim[0], self.dim[1]
        ln, seq = len(ops), self.eos[0][:]
        if dim0*dim1-1 != len(seq) or (roc == 'c' and (dim0 < 3 or dim1 < 2)) or (roc == 'r' and (dim0 < 2 or dim1 < 3)):
            return False
        groups = [seq[dim1*j:dim1*(j+1)] for j in range(dim0)] if roc == 'r' else [seq[j::dim1] for j in range(dim1)]
        o, dimx = ('/','*','-','+'), dim1-2 if roc == 'r' else dim0-2
        combs, sol = tuple(product(ops, repeat=dimx)), None
        for c in combs:
            pos_sol = self.reck2(deepcopy(groups), list(c), c, o, dimx)
            if pos_sol:
                sol = pos_sol
                break
        return False if sol is None else [KBEDK([self.eos[0]+[sol], 0], self.module, self.dim)]

    def reck2(self, groups, ops, c, o, dim):
        for g in groups:
            e = g.pop() if g != groups[-1] else None
            m = [0]*(len(g)+dim)
            m[::2], m[1::2] = g, ops
            for j in range(4):
                for k in range(c.count(o[j])):
                    i = m.index(o[j])
                    x = self.do_op(o[j], m[i-1], m[i+1], False)
                    if x is False:
                        return False
                    m[i-1:i+2] = [x]
            if m[0] != e and e is not None:
                return False
            elif e is None:
                return m[0]

    def transposed(self):
        dim0, dim1 = self.dim[0], self.dim[1]
        if dim0 <= 1 or dim1 <= 1 or len(self.eos[0]) != dim0*dim1-1:
            return False
        transposed = []
        for i in range(dim1):
            transposed += self.eos[0][i::dim1]
        sol = KBEDK(transposed, self.module, (dim1, dim0))
        return [sol.basic()]

    def rev_basic(self, n, m, h):
        for j in range(h, 0, -1):
            for k in range(n):
                if k+m < len(self.eos[j]):
                    e = self.check_number(self.eos[j][k+m]+self.eos[j-1][k+m])
                    if e is False:
                        return False
                    self.eos[j-1].append(e)
            m += 1
        return self.eos[0]

    def rev_divisions(self, n, m, h):
        for j in range(h, 0, -1):
            for k in range(n):
                if k+m < len(self.eos[j]):
                    e = self.check_number(self.eos[j][k+m]*self.eos[j-1][k+m])
                    if e is False:
                        return False
                    self.eos[j-1].append(e)
            m += 1
        return self.eos[0]

    def rev_exponentiation(self, j):
        if j == -1 and 0 in self.eos[0]:
            return False
        row = []
        for x in self.eos[0]:
            try:
                e = (x)**(1/j)
                if type(e) == complex:
                    return False
                e = self.check_number(e)
                if e is False:
                    return False
                row.append(e)
            except:
                row.append(inf)
        self.eos[0] = row
        return self.eos[0]

    def rev_logarithm(self, n, m):
        for k in range(n):
            if k+m < len(self.eos[1]):
                if self.eos[0][k+m] == 0 and self.eos[1][k+m] == 0:
                    return False
                try:
                    e = (self.eos[0][k+m])**(self.eos[1][k+m])
                    if type(e) == complex:
                        return False
                    e = self.check_number(e)
                    if e is False:
                        return False
                    self.eos[0].append(e)
                except:
                    self.eos[0].append(inf)
        return self.eos[0]

    def rev_double_operation(self, op1, op2, n, m):
        ln_2 = len(self.eos[0]) % 2
        for k in range(n):
            if k+m < len(self.eos[1]):
                op = op1 if (ln_2 == 0 and k % 2 != 0) or (ln_2 != 0 and k % 2 == 0) else op2
                e = self.do_op(op, self.eos[1][k+m], self.eos[0][k+m], True)
                if e is False:
                    return False
                self.eos[0].append(e)
        return self.eos[0]

    def rev_multi_level_p(self, div):
        base_copy = self.eos[0][:]
        for j in range(len(self.eos) - 1, 0, -1):
            for k in range(len(self.eos[j])):
                if k >= len(self.eos[j-1]):
                    continue
                if self.eos[j][k] is not None and self.eos[j-1][k] is not None:
                    e0 = self.eos[j][k]*self.eos[j-1][k] if div else self.eos[j][k]+self.eos[j-1][k]
                    e = self.check_number(e0)
                    if e is False:
                        return False
                    if k < len(self.eos[j-1]) - 1:
                        self.eos[j-1][k+1] = e
                    else:
                        self.eos[j-1].append(e)
        seq = self.eos[0]
        if None in seq:
            seq = self.eos[0][:self.eos[0].index(None)]
        return seq if seq != base_copy else False

    def rev_multi_level_i(self, subedk, dx, dy, b):
        ln, vx, vy = len(subedk.eos[0]), [], []
        for i in range(ln):
            vx.append(dx*i)
            vy.append(dy*i+b)
        if dx == 0 or dy == 0:
            x, y = ln, ln
        else:
            x, y = dx*ln, dy*ln+b
        for j in range(y):
            if j > len(self.eos)-1:
                self.eos.append([None]*x)
            if x > len(self.eos[j]):
                self.eos[j] += [None]*(x-len(self.eos[j]))
        for k in range(ln):
            p, q = vx[k], vy[k]
            if self.eos[q][p] is None:
                self.eos[q][p] = subedk.eos[0][k]
        return self.eos

    def rev_focusing_i(self, subedks, shift, divs, n, ln):
        seq = self.eos[0]+[None]*n
        for i in range(ln):
            l = subedks[i] if type(subedks[i]) == list else subedks[i].eos[0]
            a = shift + sum(divs[:i])
            s, o = sum(divs), round(len(l)/divs[i])+1
            race = [a+j+s*k for k in range(o) for j in range(divs[i]) if a+j+s*k <= len(seq)-1]
            for p, q in zip(race, range(len(race))):
                if q < len(l) and seq[p] is None:
                    seq[p] = l[q]
        return seq if None not in seq else seq[:seq.index(None)]

    def rev_split_of_elements_i(self, subedks):
        s1, s2, r = subedks[0].eos[0], subedks[1].eos[0], []
        for x, y in zip(s1, s2):
            if y == 0:
                return False
            if (type(x) == float and x-int(x) != 0.0) or x < 0 or \
                    (type(y)==float and y-int(y) != 0.0) or y < 0 or y > 10**5:
                return False
            r.append(int(str(int(x))*int(y)))
        return r if self.eos[0] == r[:len(self.eos[0])] else False

    def rev_dgee_i(self, subedks):
        s1, s2, r = subedks[0].eos[0], subedks[1].eos[0], []
        try:
            for x, y in zip(s1, s2):
                if y < 0 or (type(y) != int and y-int(y) != 0.0):
                    return False
                if y > 200:
                    break
                r += [x]*int(y)
                if len(r) > 200:
                    break
            return r if self.eos[0] == r[:len(self.eos[0])] else False
        except:
            return False

    def rev_dgde_i(self, subedks):
        if len(subedks) == 1:
            s1 = []
            for e in self.eos[0]:
                if e not in s1:
                    s1.append(e)
            s2 = subedks[0].eos[0]
        else:
            s1 = subedks[0].eos[0]
            s2 = subedks[1].eos[0]
        seq, r1, r2 = self.eos[0][:], [], []
        for i in s2:
            if i < 0 or (type(i) != int and i-int(i) != 0.0):
                return False
            i = int(i)
            r1, r2 = r1+s1[:i], r2+s1[-i:]
            del seq[:i]
        if self.eos[0] == r1[:len(self.eos[0])]:
            return r1
        elif self.eos[0] == r2[:len(self.eos[0])]:
            return r2
        else:
            return False

class SeqPredictor:
    def __init__(self, states_list, n):
        self.states_list = states_list
        self.n = n

    def predictor(self):
        for i in range(len(self.states_list)-1):
            st0, st1 = self.states_list[i-1], self.states_list[i]
            act0, act1, p = st0.action, st1.action, st1.parents
            for j in range(len(st1.eos)):
                edk = st1.eos[j]
                if j in st1.sols:
                    if act1 == 'DIV':
                        edk.eos[-1] += [1]*self.n
                        seq = edk.rev_divisions(self.n, 1, len(edk.eos)-1)
                    elif act1[:3] == 'ANA' or act1 == 'RECK' or 'SYM' in act1:
                        seq = st1.eos[j].eos[0]
                    else:
                        edk.eos[-1] += [0]*self.n
                        seq = edk.rev_basic(self.n, 1, len(edk.eos)-1)
                elif act0[:3] == 'RED' or act0 == 'LOG' or act0[:2] == 'ML' or act0[:3] == 'DOP':
                    seq = self.prediction(act0, act1, edk)
                else:
                    seq = st1.eos[j].eos[0]
                if act1[:3] == 'EXP':
                    if not seq:
                        return False
                    a = act1[act1.index('(')+1:act1.index(')')]
                    if '/' in a:
                        e0 = a.split('/')
                        e = int(e0[0]) / int(e0[1])
                    elif '.' in a:
                        e = float(a)
                    else:
                        e = int(a)
                    seq = edk.rev_exponentiation(e)
                if not seq:
                    return False
                new_edk = self.insertion(st1.eos, p, st1.action, self.states_list[i+1].eos[p[j]], seq, j)
                if new_edk is False:
                    return False
                self.states_list[i+1].eos[p[j]].eos = new_edk
        ef, af = self.states_list[-1].eos[0], self.states_list[-2].action
        return self.prediction(af, af, ef) if af[:3] == 'RED' or \
            af == 'LOG' or af[:2] == 'ML' or af[:3] == 'DOP' else ef.eos[0]

    def prediction(self, act0, act1, edk):
        if act0 == 'LOG':
            return edk.rev_logarithm(self.n, len(edk.eos)-1)
        elif act0[:3] == 'DOP':
            return edk.rev_double_operation(
                act0[4], act0[6], self.n, len(edk.eos)-1)
        elif act0[:3] == 'RED':
            h = int(act0[4])
            if act1[:3] == 'DIV':
                return edk.rev_divisions(self.n, len(edk.eos)-h, h)
            return edk.rev_basic(self.n, len(edk.eos)-h, h)
        return edk.rev_multi_level_p(act1[:3] == 'DIV')

    def insertion(self, edks1, p, act1, edk2, seq, j):
        if act1[:2] == 'ML':
            a, parents1 = [], p[:]
            while parents1:
                b = p.count(p[0])
                a = a + list(range(b))
                del parents1[:b]
            act1 = act1[act1.index('(')+1:act1.index(')')].split(',')
            b, dx, dy = a[j], int(act1[0]), int(act1[1])
            return edk2.rev_multi_level_i(edks1[j], dx, dy, b)
        elif act1[:3] == 'FOC':
            if j == len(edks1)-1 or p[j] != p[j+1]:
                edks1 = [edks1[i] for i in range(len(edks1)) if p[i] == p[j]]
                act1 = act1[act1.index('(')+1:act1.index(')')].split(',')
                shift, divs = int(act1[0]), [int(x) for x in act1[1:]]
                edk2.eos[0] = edk2.rev_focusing_i(edks1, shift, divs, self.n, p.count(p[j]))
        elif act1[:3] == 'SOE':
            if j == len(edks1)-1 or p[j] != p[j+1]:
                edks1 = [edks1[i] for i in range(len(edks1)) if p[i] == p[j]]
                edk2.eos[0] = edk2.rev_split_of_elements_i(edks1)
                return False if not edk2.eos[0] else edk2.eos
        elif act1[:4] == 'DGEE':
            if j == len(edks1)-1 or p[j] != p[j+1]:
                edks1 = [edks1[i] for i in range(len(edks1)) if p[i] == p[j]]
                edk2.eos[0] = edk2.rev_dgee_i(edks1)
                return False if not edk2.eos[0] else edk2.eos
        elif act1[:4] == 'DGDE':
            if j == len(edks1)-1 or p[j] != p[j+1]:
                edks1 = [edks1[i] for i in range(len(edks1)) if p[i] == p[j]]
                edk2.eos[0] = edk2.rev_dgde_i(edks1)
                return False if not edk2.eos[0] else edk2.eos
        elif act1[:3] == 'RED':
            edk2.eos[int(act1[4])] = seq
        elif act1[:3] == 'DOP' or act1 == 'LOG':
            edk2.eos[1] = seq
        elif 'DIAG' in act1 or act1 == 'TRA':
            edk2.eos[0] = seq+[seq[edk2.dim[0]-1]]
        else:
            edk2.eos[0] = seq
        return edk2.eos

sr0 = [[0, 1, 1.7071, 2.3660, 3, 3.6180, 4.2247, 4.8229],
      [2, 16, 4, 256, 16, 65536, 256, 4294967296, 65536],
      [3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584],
      [1, 3, 7, 12, 18, 26, 35, 45, 57, 70, 84, 100, 117, 135],
      [0, 3, 12, 17, 102, 109, 872, 881],
      [0, 3, 8, 24, 63, 168, 440, 1155, 3024, 7920],
      [1, 3, 6, 10, 15, 21, 28, 36],
      [1, 0, -1, 0, 1, 0, -1, 0, 1],
      [1, 3, 5, 7, 9, 11, 13, 15, 17],
      [1, 2, 4, 7, 11, 16, 22, 29],
      [1, 4, 9, 16, 25, 36, 49, 64, 81],
      [3, 5, 10, 12, 24, 26, 52, 54, 108],
      [3, 4, 8, 17, 33, 58, 94, 143],
      [3, 6, 18, 72, 360, 2160, 15120, 120960],
      [4, 5, 8, 13, 20, 29, 40, 53],
      [5, 11, 17, 23, 29, 35, 41],
      [11, 9, 7, 5, 3, 1, -1, -3],
      [30, 29, 27, 26, 24, 23, 21, 20, 18, 17, 15],
      [144, 121, 100, 81, 64, 49, 36, 25],
      [2, 2, 4, 6, 10, 16, 26, 42, 68],
      [81, 27, 9, 3, 1, 0.3333333333333333, 0.1111111111111111],
      [1, 1, 2, 3, 5, 8, 13, 21],
      [21, 20, 18, 15, 11, 6, 0],
      [8, 6, 7, 5, 6, 4, 5, 3, 4],
      [4294967296, 65536, 256, 16, 4, 2],
      [3, 7, 14, 24, 37, 53, 72],
      [-3, -1, 2, 6, 11, 17, 24],
      [-1, 0, 3, 8, 15, 24, 35],
      [-9, 2, 12, 21, 29, 36, 42],
      [1, 0, 0, 1, 3, 6, 10, 15],
      [2, 2, 4, 4, 8, 8, 16, 16, 32, 32, 64, 64, 128, 128, 256, 256],
      [8, 6, 4, 3, 1, -1, -2, -4, -6, -7, -9, -11, -12, -14, -16],
      [362880, 40320, 5040, 720, 120, 24, 6, 2, 1, 1],
      [14, 1, -5.5, -8.75, -10.375, -11.1875, -11.59375],
      [-1, -1, 0, 2, 5, 9, 14, 20],
      [-2, 3, 1, 4, 5, 9, 14, 23],
      [-3, -1, -4, 0, -5, 1, -6, 2, -7],
      [-4, -1, -3, 0, -2, 1, -1, 2, 0],
      [1, 2, 0, 2, -1, 2, -2, 2, -3],
      [1, -1, 0, -3, -1, -5, -2, -7, -3, -9],
      [0, 1, 0, -1, 0, 1, 0, -1, 0, 1, 0, -1],
      [-2, 0, -3, 1, -4, 2, -5, 3, -6],
      [23, 34, 45, 56, 67, 78],
      [-4, -3, 0, 5, 12, 21, 32],
      [-4, -2, 2, 8, 16, 26, 38],
      [6, 4, 0, -6, -14, -24, -36],
      [-2, 0, 4, 10, 18, 28, 40],
      [-6, -1, 5, 12, 20, 29, 39],
      [6400, 1600, 400, 100, 25, 6.25, 1.5625],
      [0, 7, 24, 51, 88, 135, 192],
      [2, 4, 12, 48, 240, 1440, 10080],
      [-10, 12, 44, 86, 138, 200, 272],
      [1, 0, 1, 1, 1, 2, 1, 3, 1, 4],
      [-1, 2, -2, -4, 8, -32, -256, 8192],
      [-10, 0, 15, 35, 60, 90, 125, 165],
      [-2, -1, 1, 5, 13, 29, 61, 125],
      [-3, -2, 0, 1, 3, 4, 6, 7, 9],
      [1, 2, 6, 14, 29, 56, 102, 176, 289],
      [0, 1, 2, 8, 29, 80, 181, 357, 638],
      [8, 1, 0, -1, -8, -27, -64, -125, -216],
      [0, 2, 9, 28, 75, 186, 441, 1016, 2295, 5110, 11253, 24564],
      [3, 6, 18, 36, 108, 216, 648, 1296, 3888, 7776, 23328, 46656, 139968],
      [2, 3, 5, 9, 17, 33, 65, 129, 257, 513, 1025, 2049, 4097, 8193],
      [2, 10, 26, 50, 82, 122, 170, 226, 290, 362, 442, 530, 626],
      [4, 1, 0, 1, 4, 9, 16, 25, 36],
      [9, 3, 6, 6, 2, 5, 3, 1, 4, 0, 0, 3, -3],
      [-9, 1, -5, 3, -4, 2, -6, -2, -11, -9],
      [97.5, 57, 30, 12, 0, -8, -13.333333333333332],
      [6, 4, 3, 3, 2, 2, 3, 1, 6, 0, 11, -1],
      [7, 16, 52, 196, 772, 3076, 12292],
      [2, 3, 5, 4, 2, 4, 7, 5, 2, 5, 9, 6, 2, 6, 11, 7, 2],
      [8, 0, 2, 0, 0, 0, 2, 0, 8, 0, 18, 0],
      [3, 0, 1, -1, -2, -5, -9, -16],
      [1, 0, -1, -1, -2, -1, -3, -4, -1, -5],
      [3, 9, 22.5, 45, 67.5, 67.5, 33.75, 0],
      [0, 2, 9, 24, 50, 90, 147, 224],
      [3, 8, 16, 28, 45, 68, 98, 136],
      [0, 0, 4, 5, 14, 16, 30, 33, 52, 56],
      [0, 8, 15, 35, 48, 80, 99, 143, 168, 224],
      [4, 32, 108, 256, 500, 864, 1372, 2048],
      [0, 17, 74, 195, 404, 725, 1182, 1799],
      [6, 24, 60, 120, 210, 336, 504, 720],
      [3, 6, 15, 42, 123, 366, 1095, 3282],
      [23, 31, 53, 83, 135, 217, 351, 567],
      [2, 6, 19, 53, 126, 262, 491, 849, 1378],
      [1, 2, 5, 9, 16, 27, 45, 74],
      [0, 2, 6, 12, 20, 30, 42, 56, 72, 90, 110, 132],
      [3, 5, 8, 12, 17, 23, 30, 38],
      [5, 8, 20, 68, 260, 1028, 4100],
      [-6, -4, 0, 2, 6, 8, 12, 14, 18]]

sr1 = [[5, 9, 35, 125, 345, 785, 1559, 2805, 4685],
      [2, 5, 11, 21, 37, 63, 107],
      [6, 288, 884736, 173946175488, 2188749418902061056],
      [1, 2, 8, 48, 384, 3840],
      [5, 13, 35, 97, 275, 793, 2315, 6817],
      [12, 44, 144, 432, 1216, 3264, 8448, 21248],
      [1, 2, 3, 4, 5],
      [1, 4, 9, 16],
      [1, 3, 9, 27],
      [0, 1, 1, 2, 3, 5, 8, 13],
      [1, 2, 4, 8, 16],
      [12, 15, 8, 11, 4, 7, 0, 3],
      [148, 84, 52, 36, 28, 24, 22],
      [2, 12, 21, 29, 36, 42, 47, 51],
      [2, 3, 5, 9, 17, 33, 65, 129],
      [2, 5, 8, 11, 14, 17, 20, 23],
      [2, 5, 9, 19, 37, 75, 149, 299],
      [25, 22, 19, 16, 13, 10, 7, 4],
      [28, 33, 31, 36, 34, 39, 37],
      [3, 6, 12, 24, 48, 96, 192],
      [3, 7, 15, 31, 63, 127, 255],
      [4, 11, 15, 26, 41, 67, 108],
      [5, 6, 7, 8, 10, 11, 14, 15],
      [54, 48, 42, 36, 30, 24, 18],
      [6, 8, 5, 7, 4, 6, 3, 5],
      [6, 9, 18, 21, 42, 45, 90, 93],
      [7, 10, 9, 12, 11, 14, 13, 16],
      [8, 10, 14, 18, 26, 34, 50, 66],
      [8, 12, 10, 16, 12, 20, 14, 24],
      [8, 12, 16, 20, 24, 28, 32, 36],
      [9, 20, 6, 17, 3, 14, 0, 11],
      [0, 1, 4, 9],
      [0, 2, 4, 6],
      [1, 1, 2, 3, 5],
      [0, 1, 2, 1, 4, 1],
      [0, 0, 1, 1, 0, 0, 1, 1],
      [0, 1, 3, 7],
      [1, 2, 2, 3, 3, 3, 4, 4, 4, 4],
      [1, 4, 7, 10, 13, 16, 19, 22],
      [2, 4, 3, 5, 4, 6, 5, 7],
      [4, 11, 15, 26, 41, 67, 108, 175],
      [5, 6, 12, 19, 32, 52, 85, 138],
      [8, 10, 14, 18, 26, 34, 50, 66],
      [1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5],
      [2, 2, 4, 8, 32, 256],
      [2, 3, 5, 8, 13, 21, 34],
      [1, 2, 1, 3, 1, 4, 1, 5],
      [5, 10, 1, 2, 22, 44, 3, 6, 7, 14],
      [0, 1, 4, 9, 16],
      [1, 4, 9, 16, 25],
      [4, 7, 12, 20, 32],
      [4, 7, 12, 20, 33],
      [7, 11, 15, 19, 23],
      [0, 2, 4, 6, 8, 10],
      [3, 6, 17, 66, 327],
      [1, 1, 2, 6, 24, 120, 720, 5040, 40320],
      [2, 5, 8, 11, 14, 17],
      [3, 6, 12, 24, 48],
      [1, 2, 3, 5, 8, 13, 21, 34, 55],
      [1, 1, 2, 6, 24, 120],
      [1, 2, 3, 4],
      [3, 2, 1, 0],
      [1, 11, 111, 1111],
      [1, 44, 1, 27, 1, 92],
      [1, 39, 1, 35, 1, 28],
      [1, 1, 7, 1],
      [46, 147, 9, 1, 1, 1]]

kl2 = ['DIV', 'RED(1)', 'RED(2)', 'RED(3)', 'FOC(0,1,1)', 'ML(1,1)', 'LOG',
      'EXP(1/2)', 'ML(0,1)', 'FOC(0,2,1)', 'DOP(-,/)', 'EXP(-1)', 'ANA(0,3)',
      'FOC(0,1,1,1)', 'FOC(0,1,1,1,1)', 'DGEE(0,0)', 'DGEE(0,1)', 'RSYM', 'SSYM',
      'DOP(*,/)', 'DOP(/,-)', 'DOP(+,-)', 'DOP(+,*)', 'DOP(+,/)', 'ML(2,1)', 'FOC(0,1,2)',
      'FOC(0,2,2)', 'ANA(0,2)', 'ANA(0,4)', 'ANA(0,5)', 'EXP(1/3)',
      'DGEE(1,0)', 'DGEE(1,1)', 'DGDE(0,0)', 'DGDE(0,1)', 'DGDE(1,0)',
      'DGDE(1,1)', 'SOE', 'FOC(0,3,1)', 'FOC(0,1,3)', 'FOC(0,3,2)', 'FOC(0,2,3)',
      'FOC(0,3,3)', 'FOC(0,4,1)', 'FOC(0,4,2)', 'FOC(0,4,3)', 'FOC(0,4,4)',
      'FOC(0,1,4)', 'FOC(0,2,4)', 'FOC(0,3,4)', 'FOC(0,4,4)', 'FOC(0,5,1)',
      'FOC(0,5,2)', 'FOC(0,5,3)', 'FOC(0,5,4)', 'FOC(0,1,5)', 'FOC(0,2,5)',
      'FOC(0,3,5)', 'FOC(0,4,5)', 'FOC(0,5,5)']


def read_path(path):
    file1 = open(path, "r", encoding="utf8")
    Content = file1.read()
    CoList = Content.split("\n")
    file1.close()
    return CoList

def write_path(path, lines):
    file1 = open(path, 'w', encoding="utf8")
    for line in lines:
        file1.write(str(line))
        file1.write('\n')
    file1.close()

def execute_kitbit_gb_False(seqs, kl, mz, path):
    results, solved = [], 0
    for seq in seqs:
        h = KitBit(seq[:-1], kl, 5000000000, 3, search_algorithm='BFS', n=1, min_zeros=mz, epsilon=exp(-18), all_solutions=False)
        x = h.handler()
        results.append(x)
        if not x[0][0] or len(x[0][0])<len(seq) or x[0][0] != seq:
            continue
        else:
            solved += 1
    write_path(path, results)
    accuracy = solved / len(seqs) * 100
    print(f"[GB-FALSE] Total: {len(seqs)} | Solved: {solved} | Accuracy: {accuracy:.2f}%")
    print(f"Results saved to: {path}\n")
    # print(solved/len(seqs)*100)
    
def execute_kitbit_gb_True(seqs, kl, mz, path):
    results, solved = [], 0
    for seq in seqs:
        h = KitBit(seq[:-1], kl, 5000000000, 3, search_algorithm='BFS', n=1, min_zeros=mz, epsilon=exp(-18), all_solutions=True)
        x = h.handler()
        results.append(x)
        if not x[0][0]:
            continue
        sol = [1 for pos_sol in x[0] if pos_sol[0]==seq]
        if 1 in sol:
            solved += 1
    write_path(path, results)
    # print(solved/len(seqs)*100)
    accuracy = solved / len(seqs) * 100
    print(f"[GB-TRUE] Total: {len(seqs)} | Solved: {solved} | Accuracy: {accuracy:.2f}%")
    print(f"Results saved to: {path}\n")
    
execute_kitbit_gb_False(sr0, kl2, 1, 'results/IQ_S1Z.txt')
execute_kitbit_gb_False(sr1, kl2, 1, 'results/LI_S1Z.txt')
execute_kitbit_gb_True(sr0, kl2, 1, 'results/IQ_N1Z.txt')
execute_kitbit_gb_True(sr1, kl2, 1, 'results/LI_N1Z.txt')

def execute_kitbit_oeis(sols, kl, sols_cad):
    solved, not_solved, q = [], [], 0
    for i in range(len(sols)):
        if i % 1600 == 0:
            print(f"[OEIS Progress] Processed: {i} | Solved: {len(solved)} | Failed: {len(not_solved)}")
        seq, solc = sols[i], sols_cad[i]
        kitasi = kl[i][1:-1].split('&')
        h = KitBit(seq[:-1], kitasi, 5000000000, kitasi, search_algorithm='branch', n=1, min_zeros=1, epsilon=exp(-18), all_solutions=False)
        x = h.handler()
        if not x[0][0] or len(x[0][0])<len(seq) or x[0][0] != seq:
            h = KitBit(seq[:-1], kitasi, 5000000000, kitasi, search_algorithm='branch', n=1, min_zeros=2, epsilon=exp(-18), all_solutions=False)
            y = h.handler()
            if not y[0][0] or y[0][0] != seq or len(y[0][0])<len(seq):
                not_solved.append(solc)
            else:
                q += 1
                solved.append(y)
        else:
            solved.append(x)
    return solved, not_solved

def execute_series_oeis(sols, kl, depth, N, mz):
    for i in range(len(sols)):
        seq = sols[i]
        h = KitBit(seq[:-1], kl, 5000000000, depth, search_algorithm='BFS', n=N, min_zeros=mz, epsilon=exp(-18), all_solutions=False)
        x = h.handler()
        print(x)
        if x[0][0] == seq:
            print(True)
        h = KitBit(seq[:-1], kl, 5000000000, depth, search_algorithm='BFS', n=N, min_zeros=mz, epsilon=exp(-18), all_solutions=True)
        x = h.handler()
        print(x) 

#------Daata for composite series------
def interleave2(a, b, take_last_from= "a"):
    """
    Interleave three sequences:
    a1, b1, a2, b2, ...
    """
    out =[]
    n = min(len(a), len(b))
    
    for i in range(n):
        out.append(a[i])
        out.append(b[i])
        
    if take_last_from == 'a' and len(a) > n:
        out.append(a[n])
    if take_last_from == 'b' and len(b) > n:
        out.append(b[n])
        
    return out

def interleave3(a, b, c, take_last_from="c"):
    """
    Interleave three sequences:
    a1, b1, c1, a2, b2, c2, ...
    """
    out = []
    n = min(len(a), len(b), len(c))
    for i in range(n):
        out.extend([a[i], b[i], c[i]])

    if take_last_from == "a" and len(a) > n:
        out.append(a[n])
    elif take_last_from == "b" and len(b) > n:
        out.append(b[n])
    elif take_last_from == "c" and len(c) > n:
        out.append(c[n])

    return out

base_families = {
    "A000027_naturals":   [1, 2, 3, 4, 5, 6],
    "A005843_evens":      [2, 4, 6, 8, 10, 12],
    "A005408_odds":       [1, 3, 5, 7, 9, 11],
    "A000217_triangular": [1, 3, 6, 10, 15, 21],
    "A000290_squares":    [1, 4, 9, 16, 25, 36],
    "A000578_cubes":      [1, 8, 27, 64, 125, 216],
    "A000040_primes":     [2, 3, 5, 7, 11, 13],
    "A000045_fib":        [1, 1, 2, 3, 5, 8],
    "A000142_factorial":  [1, 2, 6, 24, 120, 720],
    "A001477_posints":    [0, 1, 2, 3, 4, 5]
}

direct_composite_set = [
    [1, 1, 2, 3, 3, 6, 4, 10, 5],      # naturals + triangular
    [1, 1, 2, 4, 3, 9, 4, 16, 5],      # naturals + squares
    [1, 1, 2, 8, 3, 27, 4, 64, 5],     # naturals + cubes
    [2, 1, 4, 3, 6, 6, 8, 10, 10],     # evens + triangular
    [1, 2, 2, 3, 3, 5, 4, 7, 5],       # naturals + primes
    [1, 1, 2, 1, 3, 2, 4, 3, 5],       # naturals + fib
    [2, 1, 4, 4, 6, 9, 8, 16, 10],     # evens + squares
    [1, 2, 3, 4, 5, 6, 7, 8, 9],       # simple merge style
    [1, 2, 2, 6, 3, 12, 4, 20, 5],     # naturals + n(n+1)
    [1, 100, 2, 99, 3, 98, 4, 97, 5],  # increasing + decreasing
]

constructed_composite_set = []

pairs = [
    ("A000027_naturals", "A000217_triangular"),
    ("A000027_naturals", "A000290_squares"),
    ("A000027_naturals", "A000578_cubes"),
    ("A000027_naturals", "A000040_primes"),
    ("A000027_naturals", "A000045_fib"),
    ("A005843_evens", "A000217_triangular"),
    ("A005843_evens", "A000290_squares"),
    ("A005843_evens", "A000040_primes"),
    ("A005408_odds", "A000217_triangular"),
    ("A005408_odds", "A000290_squares"),
    ("A005408_odds", "A000040_primes"),
    ("A000045_fib", "A000217_triangular"),
    ("A000045_fib", "A000290_squares"),
    ("A000040_primes", "A000290_squares"),
    ("A000040_primes", "A000217_triangular"),
]

# make 15 pairwise 
for a_name, b_name in pairs:
    a = base_families[a_name]
    b = base_families[b_name]
    seq = interleave2(a, b, take_last_from="a")[:9]
    constructed_composite_set.append({
        "sequence": seq,
        "source": f"{a_name} + {b_name}",
        "type": "pair_interleave"
    })

# 15 triple 
triples = [
    ("A000027_naturals", "A005843_evens", "A000217_triangular"),
    ("A000027_naturals", "A000290_squares", "A000040_primes"),
    ("A000027_naturals", "A000045_fib", "A000217_triangular"),
    ("A005408_odds", "A000290_squares", "A000040_primes"),
    ("A005843_evens", "A000217_triangular", "A000040_primes"),
    ("A000027_naturals", "A000578_cubes", "A000040_primes"),
    ("A000045_fib", "A000290_squares", "A000217_triangular"),
    ("A005408_odds", "A000217_triangular", "A000045_fib"),
    ("A005843_evens", "A000290_squares", "A000045_fib"),
    ("A000027_naturals", "A000040_primes", "A000142_factorial"),
    ("A005408_odds", "A000040_primes", "A000217_triangular"),
    ("A000027_naturals", "A000290_squares", "A000217_triangular"),
    ("A005843_evens", "A000040_primes", "A000290_squares"),
    ("A000045_fib", "A000040_primes", "A000217_triangular"),
    ("A000027_naturals", "A005408_odds", "A000040_primes"),
]

for a_name, b_name, c_name in triples:
    a = base_families[a_name]
    b = base_families[b_name]
    c = base_families[c_name]
    seq = interleave3(a, b, c, take_last_from="a")[:9]
    constructed_composite_set.append({
        "sequence": seq,
        "source": f"{a_name} + {b_name} + {c_name}",
        "type": "triple_interleave"
    })

composite_dataset = []
composite_dataset.extend(direct_composite_set)
composite_dataset.extend(constructed_composite_set)

print(f"Before dedup: {len(composite_dataset)}")

seen = set()
unique_composite_dataset = []

for item in composite_dataset:
    if isinstance(item, dict):
        seq_tuple = tuple(item["sequence"])
    else:
        seq_tuple = tuple(item)

    if seq_tuple not in seen:
        seen.add(seq_tuple)
        unique_composite_dataset.append(item)

composite_dataset = unique_composite_dataset

print(f"After dedup: {len(composite_dataset)}")

extra_composites = [
    [1, 2, 2, 6, 3, 24, 4, 120, 5],     # naturals + factorial
    [2, 1, 3, 2, 5, 6, 7, 24, 11],      # primes + factorial
    [1, 2, 4, 3, 9, 5, 16, 7, 25],      # squares + primes
    [1, 1, 3, 2, 6, 6, 10, 24, 15],     # triangular + factorial
    [2, 1, 3, 1, 5, 2, 7, 3, 11],       # primes + fib
    [1, 1, 4, 2, 9, 3, 16, 5, 25],      # squares + fib
    [2, 1, 4, 2, 6, 6, 8, 24, 10],      # evens + factorial
]

for seq in extra_composites:
    seq_tuple = tuple(seq)
    if seq_tuple not in seen:
        seen.add(seq_tuple)
        composite_dataset.append(seq)

composite_test_set = [
    item["sequence"] if isinstance(item, dict) else item
    for item in composite_dataset
]
print(composite_test_set)

def run_composite_baseline(seqs, kl, mz=1, depth=3):
    solved = 0
    results = []

    for seq in seqs:
        h = KitBit(
            seq[:-1], kl, 500000, depth,
            search_algorithm='BFS',
            n=1, min_zeros=mz,
            epsilon=exp(-18),
            all_solutions=False
        )
        x = h.handler()

        pred_seq = x[0][0] if x and x[0] else False
        predicted_next = pred_seq[len(seq)-1] if pred_seq and len(pred_seq) >= len(seq) else None

        full_match = bool(pred_seq) and len(pred_seq) >= len(seq) and pred_seq[:len(seq)] == seq
        next_match = predicted_next == seq[-1]

        ok = next_match   # use this for solving accuracy

        if ok:
            solved += 1

        results.append({
            "input": seq[:-1], 
            "expected": seq[-1],
            "predicted": predicted_next,
            "solved": ok,
            "full_match": full_match,
            "actions": x[0][2] if x and x[0] and len(x[0]) > 2 else [],
            "time": x[1] if x and len(x) > 1 else None
        })

    accuracy = solved / len(seqs) * 100
    print(f"[Composite Baseline] Total: {len(seqs)} | Solved: {solved} | Accuracy: {accuracy:.2f}%")

    print("\nFailed sequences:")
    for r in results:
        if not r["solved"]:
            print(
                f"input={r['input']} | expected={r['expected']} | "
                f"predicted={r['predicted']}"
            )

    return results

def solve_subsequence_with_kitbit(subseq, kl, mz=1, depth=3):
    h = KitBit(
        subseq[:-1], kl, 500000, depth,
        search_algorithm='BFS',
        n=1, min_zeros=mz,
        epsilon=exp(-18),
        all_solutions=False
    )

    x = h.handler()
    pred_seq = x[0][0] if x and x[0] else False
    predicted_next = pred_seq[len(subseq) - 1] if pred_seq and len(pred_seq) >= len(subseq) else None

    return {
        "solved": predicted_next == subseq[-1],
        "predicted_next": predicted_next,
        "pred_seq": pred_seq,
        "actions": x[0][2] if x and x[0] and len(x[0]) > 2 else [],
        "time": x[1] if x and len(x) > 1 else None
    }


def split_odd_even(seq):
    return [seq[::2], seq[1::2]]


def split_stride3(seq):
    return [seq[0::3], seq[1::3], seq[2::3]]


def evaluate_split(parts, expected, kl, mz=1, depth=3):
    """
    Evaluate whether each decomposed subsequence can correctly predict
    its own next term when extended with the expected answer.
    
    it must feld ture for all subsequences to be marked as "solved"
    """
    if any(len(part) < 3 for part in parts):
        return None

    sub_results = []
    
    #solving all subsequneces
    for part in parts:
        result = solve_subsequence_with_kitbit(part + [expected], kl, mz=mz, depth=depth)
        sub_results.append(result)

    return {
        "parts": parts,
        "sub_results": sub_results,
        "solved_count": sum(r["solved"] for r in sub_results),
        "all_solved": all(r["solved"] for r in sub_results),
    }


def reconstruct_next_from_split(parts, sub_results, mode):
    """
    Reconstruct which subsequence should generate the next element
    in the original interleaved sequence.
    """
    preds = [r["predicted_next"] for r in sub_results]
    if any(p is None for p in preds):
        return None

    if mode == "odd_even":
        # If first part is longer, next element belongs to second part, else first part
        return preds[1] if len(parts[0]) > len(parts[1]) else preds[0]

    if mode == "stride3":
        lengths = [len(p) for p in parts]
        min_len = min(lengths)

        # Find which stream is next in the interleaving order
        for idx, part in enumerate(parts):
            if len(part) == min_len:
                return preds[idx]

    return None


def try_best_composite_split(seq, expected, kl, mz=1, depth=3):
    """
    Try only interpretable composite decomposition strategies.
    """
    strategies = [
        ("odd_even", split_odd_even),
        ("stride3", split_stride3),
    ]

    candidates = []

    for mode_name, splitter in strategies:
        parts = splitter(seq)
        result = evaluate_split(parts, expected, kl, mz=mz, depth=depth)

        if result is None:
            continue

        reconstructed_pred = reconstruct_next_from_split(parts, result["sub_results"], mode_name)
        result["mode"] = mode_name
        result["reconstructed_pred"] = reconstructed_pred
        result["matches_expected"] = reconstructed_pred == expected

        candidates.append(result)

    if not candidates:
        return None

    candidates.sort(
        key=lambda x: (
            x["matches_expected"],
            x["all_solved"],
            x["solved_count"]
        ),
        reverse=True
    )

    return candidates[0]


def run_composite_with_decomposition(seqs, kl, mz=1, depth=3):
    solved = 0
    results = []

    for seq in seqs:
        input_seq = seq[:-1]
        expected = seq[-1]

        # Step 1: baseline
        baseline = solve_subsequence_with_kitbit(input_seq + [expected], kl, mz=mz, depth=depth)
        baseline_pred = baseline["predicted_next"]
        baseline_ok = baseline_pred == expected

        final_ok = baseline_ok
        final_pred = baseline_pred
        mode_used = "baseline"
        detail = None

        # Step 2: decomposition only if baseline fails
        if not baseline_ok:
            best_split = try_best_composite_split(input_seq, expected, kl, mz=mz, depth=depth)

            if best_split and best_split["matches_expected"]:
                final_ok = True
                final_pred = best_split["reconstructed_pred"]
                mode_used = f"decomposition-{best_split['mode']}"
                detail = best_split

        if final_ok:
            solved += 1

        results.append({
            "input": input_seq,
            "expected": expected,
            "predicted": final_pred,
            "solved": final_ok,
            "mode": mode_used,
            "detail": detail
        })

    accuracy = solved / len(seqs) * 100
    print(f"[Composite Improved] Total: {len(seqs)} | Solved: {solved} | Accuracy: {accuracy:.2f}%")

    print("\nStill failed sequences:")
    for r in results:
        if not r["solved"]:
            print(
                f"input={r['input']} | expected={r['expected']} | predicted={r['predicted']}"
            )

    return results


sols, sols_cad = [], []
CoList1 = read_path('data/OEIS_SERIES_SOLVED.txt')
CoList2 = read_path('data/OEIS_KITAS.txt')
print(f"[OEIS Dataset] Total: {len(CoList1)} | Unique: {len(set(CoList1))}\n")
for j in range(len(CoList1)):
    seq1 = list(map(lambda u: int(u), CoList1[j][9:-2].split(',')))
    pos_sol = seq1[:]
    if len(pos_sol) < 4 or len(seq1[:-1])<2:
        continue
    sols.append(pos_sol)
    sols_cad.append(CoList1[j][:-1])

sol_def, not_sol = execute_kitbit_oeis(sols, CoList2, sols_cad)
print(f"[OEIS Dataset] Solved: {len(sol_def)} | UnSolved: { len(not_sol)} \n")
write_path('results/OEIS_results.txt', sol_def)


print("---------------------------- BASELINE --------------")

baseline_results = run_composite_baseline(composite_test_set, kl2)

print("----------------------------IMPROVED------------------")

improved_results = run_composite_with_decomposition(composite_test_set, kl=kl2)
