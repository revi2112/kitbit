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


class ParameterTuner:
    """Tunes KitBit algorithm parameters for optimal performance"""
    
    def __init__(self, sequences, kl, mz=1, module=10):
        self.sequences = sequences
        self.kl = kl
        self.mz = mz
        self.module = module
    
    def tune_depth(self, depths=[2, 3, 4, 5], sample_size=None):
        """Test different tree depths and return performance metrics
        
        Args:
            depths: List of depths to test
            sample_size: Number of sequences to sample (None = all)
        
        Returns:
            Dictionary: {depth: {accuracy, avg_time, solved}}
        """
        results = {}
        test_seqs = self.sequences[:sample_size] if sample_size else self.sequences
        
        for depth in depths:
            solved = 0
            times = []
            
            for seq in test_seqs:
                if len(seq) < 2:
                    continue
                
                kb = KitBit(seq, self.kl, self.mz, depth, module=self.module)
                start = default_timer()
                result = kb.handler()
                end = default_timer()
                times.append(end - start)
                
                if result[0] is not False and result[0] is not None:
                    solved += 1
            
            accuracy = (solved / len(test_seqs) * 100) if test_seqs else 0
            avg_time = sum(times) / len(times) if times else 0
            
            results[depth] = {
                'accuracy': accuracy,
                'avg_time': avg_time,
                'solved': solved
            }
        
        return results
    
    def tune_epsilon(self, epsilons=[exp(-18), exp(-15), exp(-12)], depth=3, sample_size=None):
        """Test different epsilon values for convergence threshold
        
        Args:
            epsilons: List of epsilon values to test
            depth: Tree depth to use
            sample_size: Number of sequences to sample
        
        Returns:
            Dictionary: {epsilon: {accuracy, solved}}
        """
        results = {}
        test_seqs = self.sequences[:sample_size] if sample_size else self.sequences
        
        for eps in epsilons:
            solved = 0
            
            for seq in test_seqs:
                if len(seq) < 2:
                    continue
                
                kb = KitBit(seq, self.kl, self.mz, depth, epsilon=eps, module=self.module)
                result = kb.handler()
                
                if result[0] is not False and result[0] is not None:
                    solved += 1
            
            accuracy = (solved / len(test_seqs) * 100) if test_seqs else 0
            results[eps] = {
                'accuracy': accuracy,
                'solved': solved
            }
        
        return results


class KitBitWithEarlyStopping(KitBit):
    """Extended KitBit with early stopping for improved efficiency"""
    
    def __init__(self, structure, kl, mni, depth, search_algorithm='BFS', n=2, 
                 min_zeros=1, epsilon=exp(-18), all_solutions=False,
                 max_iterations=None, early_stop_threshold=0.8):
        super().__init__(structure, kl, mni, depth, search_algorithm, n, min_zeros, epsilon, all_solutions)
        self.max_iterations = max_iterations
        self.early_stop_threshold = early_stop_threshold
        self.iteration_count = 0
    
    def numeric_solver_early_stop(self, seq, kl, depth, module):
        """Single solution solver with early stopping"""
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
        alg.max_iterations = self.max_iterations
        alg.early_stop_threshold = self.early_stop_threshold
        
        if self.search_alg == 'BFS':
            lst = alg.bfs_early_stop()
        elif self.search_alg == 'DFS':
            lst = alg.dfs([st0])
        else:
            lst = alg.branch()
        
        if lst is None or lst is False or len(lst[0]) == 1:
            return False, 0, lst
        return lst[0][0], lst[1], lst[0][1]
    
    def numeric_solver_all_sols_early_stop(self, seq, kl, depth, module):
        """All solutions solver with early stopping"""
        edk0 = KBEDK(seq[:], module, None)
        edk0 = edk0.basic()
        if edk0.is_goal('BASIC', self.epsilon, self.min_zeros, 0):
            edk0.eos[-1] += [0]*self.n
            return edk0.rev_basic(self.n, 1, len(edk0.eos)-1), 0, ['BASIC']
        
        st0 = SeqState(None, [edk0], [], [])
        alg = SeqSearchAlgorithm(st0, kl, self.mni, depth, self.n, self.min_zeros, self.epsilon)
        alg.max_iterations = self.max_iterations
        alg.early_stop_threshold = self.early_stop_threshold
        
        if self.search_alg == 'BFS':
            lst = alg.bfs_early_stop()
        elif self.search_alg == 'DFS':
            lst = alg.dfs([st0])
        else:
            lst = alg.branch()
        
        if lst is None or lst is False or len(lst[0]) == 1:
            return False, 0, lst
        return lst[0][0], lst[1], lst[0][1]


class SeqSearchAlgorithm_EarlyStop(SeqSearchAlgorithm):
    """Extension of SeqSearchAlgorithm with early stopping support"""
    
    def __init__(self, initial_state, kl, mni, depth, n=2, min_zeros=1, epsilon=exp(-18)):
        super().__init__(initial_state, kl, mni, depth, n, min_zeros, epsilon)
        self.max_iterations = None
        self.early_stop_threshold = 0.8
        self.iteration_count = 0
    
    def bfs_early_stop(self):
        """BFS with early stopping capability"""
        if self.initial_state is None:
            return None
        
        visited_states = set()
        queue = [self.initial_state]
        iteration = 0
        
        while queue and (self.max_iterations is None or iteration < self.max_iterations):
            iteration += 1
            new_queue = []
            
            for state in queue:
                state_id = id(state)
                if state_id in visited_states:
                    continue
                visited_states.add(state_id)
                
                # Check if goal is reached
                if state.is_goal():
                    return state.get_solution()
                
                # Generate successors
                successors = state.expand(self.kl, self.mni, self.depth - 1, self.n, self.min_zeros, self.epsilon, self.module)
                
                if successors:
                    for successor in successors:
                        if id(successor) not in visited_states:
                            new_queue.append(successor)
            
            queue = new_queue
        
        return None
