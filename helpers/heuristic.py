from math import exp
from collections import defaultdict

from helpers.generic import KBEDK, DynamicKitBit, KitBit, SeqSearchAlgorithm, SeqState


def learn_kita_success_rates(seqs, kl, mz=1, depth=3, epsilon=exp(-18), n=1):
    success_counts = defaultdict(int)
    for seq in seqs:
        h = KitBit(seq[:-1], kl, 5_000_000, depth,
                   search_algorithm='BFS', n=n,
                   min_zeros=mz, epsilon=epsilon, all_solutions=False)
        result = h.handler()
        if not result or not result[0] or not result[0][0]:
            continue
        pred_seq, _, actions = result[0]
        if len(pred_seq) < len(seq) or pred_seq != seq:
            continue
        for act in actions:
            if act and act != 'BASIC':
                success_counts[act] += 1
    return success_counts


def make_heuristic_kl(kl, success_counts):
    # list of seen and unseen kita 
    seen = sorted([k for k in  kl if success_counts.get(k,0) > 0], key= lambda k: -success_counts[k] )
    unseen = [k for k in kl if success_counts.get(k, 0) == 0]
    return seen + unseen


class InstrumentedBFS:
    def __init__(self, init_state, kl, mni, depth, n, min_zeros, epsilon):
        self.init_state = init_state
        self.kl = kl
        self.mni = mni
        self.depth = depth
        self.n = n
        self.min_zeros = min_zeros
        self.epsilon = epsilon

    def run(self):
        road = [[self.init_state]]
        i, k, ln = 0, -1, len(self.kl) - 1
        n_iter = [len(self.kl) ** p for p in range(1, self.depth + 1)]
        nodes = 0
        for j in range(1, sum(n_iter) + 1):
            k = k + 1 if k < ln else 0
            if not road[i]:
                st = False
            else:
                nodes += 1
                try:
                    st = road[i][-1].new_state(self.n, road[i][-1].action,
                                               self.kl[k], self.epsilon, self.min_zeros)
                except Exception:
                    st = False
            if st is False:
                i = i + 1 if k == ln else i
                road.append(False)
                continue
            road.append(road[i] + [st])
            i = i + 1 if k == ln else i
            if len(st.sols) == len(st.eos):
                return road[-1], j, nodes
            if j >= self.mni:
                return False, j, nodes
        return False, sum(n_iter), nodes


def count_nodes(seq, kl, mz, depth, epsilon=exp(-18), n=1):
    edk0 = KBEDK(seq[:-1], 10, None).basic()
    if edk0 is False:
        return 0
    if edk0.is_goal('BASIC', epsilon, mz, 0):
        return 1
    st0 = SeqState(None, [edk0], [], [])
    bfs = InstrumentedBFS(st0, kl, 5_000_000, depth, n, mz, epsilon)
    _, _, nodes = bfs.run()
    return nodes


def solve(seq, kl, mz, depth, epsilon=exp(-18), n=1):
    h = KitBit(seq[:-1], kl, 5000000, depth,
               search_algorithm='BFS', n=n,
               min_zeros=mz, epsilon=epsilon, all_solutions=False)
    r = h.handler()
    return (r and r[0] and r[0][0]
            and len(r[0][0]) >= len(seq)
            and r[0][0] == seq)


def run_heuristic_search(sr0, sr1, kl2):

    # --- Phase 1: learn ---
    print("\nPhase 1: learning kita success rates from sr0\n")
    counts = learn_kita_success_rates(sr0, kl2, mz=1, depth=3)
    kl_heur = make_heuristic_kl(kl2, counts)
    
    # sr0_train, sr0_test = sr0[:len(sr0)//2], sr0[len(sr0)//2:]
    # counts   = learn_kita_success_rates(sr0_train, kl2, ...)
    # kl_heur  = make_heuristic_kl(kl2, counts)

    for kita, cnt in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {kita:<30} {cnt} hits  ({cnt/len(sr0)*100:.1f}%)")

    print(f"\nheuristic order (top 10): {kl_heur[:10]}\n")

    # --- Phase 2: compare ---
    print("Phase 2: before vs after\n")

    for label, seqs in [("sr0 - IQ series (train, seen)", sr0), ("sr1 - literature (benchmark, un seen)", sr1)]:
        solved_before = solved_after = 0
        nodes_before_list = []
        nodes_after_list  = []
        faster = slower = same = fallbacks = 0

        for seq in seqs:
            nb = count_nodes(seq, kl2, mz=1, depth=3)
            na = count_nodes(seq, kl_heur, mz=1, depth=3)

            sb = solve(seq, kl2, mz=1, depth=3)
            sa = solve(seq, kl_heur, mz=1, depth=3)

            if sb and not sa:
                # heuristic failed, fell back to baseline — count it as baseline cost
                fallbacks += 1

            solved_before += sb
            solved_after  += sa
            nodes_before_list.append(nb)
            nodes_after_list.append(na)

            if na < nb:   faster += 1
            elif na > nb: slower += 1
            else:         same   += 1

        avg_nb = sum(nodes_before_list) / len(seqs)
        avg_na = sum(nodes_after_list)  / len(seqs)
        saved  = sum(b - a for b, a in zip(nodes_before_list, nodes_after_list) if a < b)

        print(f"{label}  ({len(seqs)} seqs)")
        print(f"  accuracy:      {solved_before}/{len(seqs)} -> {solved_after}/{len(seqs)}")
        print(f"  avg nodes:     {avg_nb:.1f} -> {avg_na:.1f}  ({(avg_na-avg_nb)/avg_nb*100:+.1f}%)")
        print(f"  faster:        {faster} seqs  ({saved} nodes saved)")
        print(f"  slower:        {slower} seqs")
        print(f"  unchanged:     {same} seqs")
        print(f"  regressions:   {fallbacks}")
        print()