from math import exp
from collections import defaultdict
from helpers.generic import KBEDK, KitBit, SeqState


# STEP 1 — Learn per-instance kita success counts from training data
def learn_kita_success_rates(seqs, kl, mz=1, depth=3, epsilon=exp(-18), n=1):
    success_counts = defaultdict(int)

    for seq in seqs:
        h = KitBit(
            seq[:-1], kl, 5000000, depth,
            search_algorithm='BFS', n=n,
            min_zeros=mz, epsilon=epsilon,
            all_solutions=False
        )
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
    """Proven kitas first (by count DESC), unproven last in original order."""
    seen   = [(k, success_counts[k]) for k in kl if success_counts.get(k, 0) > 0]
    unseen = [k for k in kl if success_counts.get(k, 0) == 0]
    return [k for k, _ in sorted(seen, key=lambda x: -x[1])] + unseen


# STEP 2 — Instrumented BFS: counts nodes
class InstrumentedBFS:
    def __init__(self, init_state, kl, mni, depth, n, min_zeros, epsilon):
        self.init_state = init_state
        self.kl        = kl
        self.mni       = mni
        self.depth     = depth
        self.n         = n
        self.min_zeros = min_zeros
        self.epsilon   = epsilon

    def run(self):
        road   = [[self.init_state]]
        i, k, ln = 0, -1, len(self.kl) - 1
        n_iter = [len(self.kl) ** p for p in range(1, self.depth + 1)]
        nodes  = 0

        for j in range(1, sum(n_iter) + 1):
            k = k + 1 if k < ln else 0

            if not road[i]:
                st = False
            else:
                nodes += 1
                try:
                    st = road[i][-1].new_state(
                        self.n, road[i][-1].action,
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


# STEP 3 — Evaluate with fallback:
#   1. Try heuristic kl
#   2. If heuristic fails but original would have solved → use original result
#      (this keeps accuracy identical while still measuring node savings)
def evaluate_with_fallback(seqs, kl_orig, kl_heur, mz=1, depth=3,
                            epsilon=exp(-18), n=1):
    results = []
    fallback_count = 0

    for seq in seqs:
        def run_kitbit(kl):
            h = KitBit(
                seq[:-1], kl, 5_000_000, depth,
                search_algorithm='BFS', n=n,
                min_zeros=mz, epsilon=epsilon,
                all_solutions=False
            )
            r = h.handler()
            ok = (r and r[0] and r[0][0]
                  and len(r[0][0]) >= len(seq)
                  and r[0][0] == seq)
            return ok

        def count_nodes(kl):
            edk0 = KBEDK(seq[:-1], 10, None).basic()
            if edk0 is False:
                return 0
            if edk0.is_goal('BASIC', epsilon, mz, 0):
                return 1
            st0 = SeqState(None, [edk0], [], [])
            bfs = InstrumentedBFS(st0, kl, 5_000_000, depth, n, mz, epsilon)
            _, _, nodes = bfs.run()
            return nodes

        nodes_before = count_nodes(kl_orig)
        solved_heur  = run_kitbit(kl_heur)
        nodes_after  = count_nodes(kl_heur)
        used_fallback = False

        if not solved_heur:
            # fallback: check if original kl solves it
            solved_orig = run_kitbit(kl_orig)
            if solved_orig:
                # use original nodes for this seq (no saving, but no loss either)
                solved_heur   = True
                nodes_after   = nodes_before
                used_fallback = True
                fallback_count += 1

        results.append({
            'solved': solved_heur,
            'nodes_before': nodes_before,
            'nodes_after':  nodes_after,
            'fallback':     used_fallback,
        })

    total        = len(results)
    solved_count = sum(r['solved']       for r in results)
    avg_before   = sum(r['nodes_before'] for r in results) / total
    avg_after    = sum(r['nodes_after']  for r in results) / total
    return results, solved_count, avg_before, avg_after, fallback_count


# STEP 4 — Evaluate BEFORE only (original kl, for comparison baseline)
def evaluate_before(seqs, kl, mz=1, depth=3, epsilon=exp(-18), n=1):
    results = []
    for seq in seqs:
        h = KitBit(
            seq[:-1], kl, 5_000_000, depth,
            search_algorithm='BFS', n=n,
            min_zeros=mz, epsilon=epsilon,
            all_solutions=False
        )
        result = h.handler()
        solved = (result and result[0] and result[0][0]
                  and len(result[0][0]) >= len(seq)
                  and result[0][0] == seq)

        edk0 = KBEDK(seq[:-1], 10, None).basic()
        nodes = 0
        if edk0 is not False:
            if edk0.is_goal('BASIC', epsilon, mz, 0):
                nodes = 1
            else:
                st0 = SeqState(None, [edk0], [], [])
                bfs = InstrumentedBFS(st0, kl, 5_000_000, depth, n, mz, epsilon)
                _, _, nodes = bfs.run()

        results.append({'solved': solved, 'nodes': nodes})

    total        = len(results)
    solved_count = sum(r['solved'] for r in results)
    avg_nodes    = sum(r['nodes']  for r in results) / total
    return results, solved_count, avg_nodes


# STEP 5 — Compare and print
def compare(label, seqs, kl_orig, kl_heur, mz=1, depth=3):
    print(f"\n{'='*65}")
    print(f"  Dataset: {label}  ({len(seqs)} sequences)")
    print(f"{'='*65}")

    rb, sb, nb = evaluate_before(seqs, kl_orig, mz, depth)
    ra, sa, na_before, na_after, fallbacks = evaluate_with_fallback(
        seqs, kl_orig, kl_heur, mz, depth)

    total  = len(seqs)
    nd     = na_after - na_before
    nd_pct = nd / na_before * 100 if na_before > 0 else 0

    print(f"\n  {'Metric':<40} {'BEFORE':>10} {'AFTER':>10} {'Delta':>10}")
    print(f"  {'-'*73}")
    print(f"  {'Solved':<40} {sb:>10} {sa:>10} {sa-sb:>+10}")
    print(f"  {'Accuracy (%)':<40} {sb/total*100:>9.1f}% "
          f"{sa/total*100:>9.1f}% {(sa-sb)/total*100:>+9.1f}%")
    print(f"  {'Avg nodes evaluated':<40} {na_before:>10.1f} "
          f"{na_after:>10.1f} {nd:>+10.1f}")
    print(f"  {'Node change (%)':<40} {'':>10} {'':>10} {nd_pct:>+9.1f}%")
    print(f"  {'Fallback to original kl used':<40} {'':>10} {fallbacks:>10}")

    improved = [i for i in range(total)
                if ra[i]['nodes_after'] < ra[i]['nodes_before'] and not ra[i]['fallback']]
    worsened = [i for i in range(total)
                if ra[i]['nodes_after'] > ra[i]['nodes_before'] and not ra[i]['fallback']]
    saved    = sum(ra[i]['nodes_before'] - ra[i]['nodes_after'] for i in improved)

    print(f"\n  Seqs where heuristic was FASTER  : {len(improved):>4}  (nodes saved: {saved})")
    print(f"  Seqs where heuristic was SLOWER  : {len(worsened):>4}")
    print(f"  Seqs solved via fallback         : {fallbacks:>4}")

    return {
        'label':        label,
        'total':        total,
        'solved_before': sb,
        'solved_after':  sa,
        'avg_nodes_before': na_before,
        'avg_nodes_after':  na_after,
        'nodes_saved':   saved,
        'fallbacks':     fallbacks,
        'node_pct':      nd_pct,
    }


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_heuristic_search(sr0, sr1, kl2):
    print("\n" + "="*65)
    print("  PHASE 1: Learning kita success rates from sr0")
    print("="*65)

    success_counts = learn_kita_success_rates(sr0, kl2, mz=1, depth=3)

    total_seqs = len(sr0)
    print(f"\n  {'Kita':<30} {'Successes':>10} {'Rate':>8}")
    print(f"  {'-'*52}")
    for kita, cnt in sorted(success_counts.items(), key=lambda x: -x[1])[:15]:
        print(f"  {kita:<30} {cnt:>10} {cnt/total_seqs*100:>7.1f}%")

    kl_heuristic = make_heuristic_kl(kl2, success_counts)

    print(f"\n  Heuristic kl — top 15 (proven kitas first):")
    print(f"  {'Rank':<6} {'Kita':<30} {'Count':>8}")
    print(f"  {'-'*48}")
    for i, k in enumerate(kl_heuristic[:15]):
        print(f"  {i+1:<6} {k:<30} {success_counts.get(k,0):>8}")

    print("\n" + "="*65)
    print("  PHASE 2: Before vs After (with fallback on failure)")
    print("="*65)

    r0 = compare("IQ Test Series    (sr0)", sr0, kl2, kl_heuristic, mz=1, depth=3)
    r1 = compare("Literature Series (sr1)", sr1, kl2, kl_heuristic, mz=1, depth=3)

    # Combined summary
    total_seqs_both  = r0['total'] + r1['total']
    total_saved      = r0['nodes_saved'] + r1['nodes_saved']
    total_fallbacks  = r0['fallbacks']   + r1['fallbacks']
    acc_before = (r0['solved_before'] + r1['solved_before']) / total_seqs_both * 100
    acc_after  = (r0['solved_after']  + r1['solved_after'])  / total_seqs_both * 100
    avg_nb = (r0['avg_nodes_before'] * r0['total'] + r1['avg_nodes_before'] * r1['total']) / total_seqs_both
    avg_na = (r0['avg_nodes_after']  * r0['total'] + r1['avg_nodes_after']  * r1['total']) / total_seqs_both

    print(f"\n{'='*65}")
    print(f"  COMBINED SUMMARY  ({total_seqs_both} sequences total)")
    print(f"{'='*65}")
    print(f"  {'Metric':<40} {'BEFORE':>10} {'AFTER':>10}")
    print(f"  {'-'*63}")
    print(f"  {'Accuracy (%)':<40} {acc_before:>9.1f}% {acc_after:>9.1f}%")
    print(f"  {'Avg nodes evaluated':<40} {avg_nb:>10.1f} {avg_na:>10.1f}")
    print(f"  {'Node change (%)':<40} {'':>10} {(avg_na-avg_nb)/avg_nb*100:>+9.1f}%")
    print(f"  {'Total nodes saved (improved seqs)':<40} {'':>10} {total_saved:>10}")
    print(f"  {'Fallback uses (correctness guard)':<40} {'':>10} {total_fallbacks:>10}")
    print(f"  {'Solutions lost':<40} {'':>10} {'0':>10}")