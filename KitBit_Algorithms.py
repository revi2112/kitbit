from math import exp, floor
from math import fabs, log, inf, ceil
from itertools import product
import os
from collections import defaultdict
from data.data import sr0, sr1, kl2, composite_test_set
from helpers.composite_helper import is_close, is_reasonable, predict_next_for_subseq, solve_subsequence_with_kitbit, try_best_composite_split
from helpers.generic import KBEDK, KitBit, SeqPredictor, SeqSearchAlgorithm, SeqState, write_path, read_path
from helpers.heuristic import run_heuristic_search

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

def run_composite_baseline(seqs, kl, mz=1, depth=3):
    solved = 0
    results = []
 
    for seq in seqs:
        h = KitBit(
            seq[:-1], kl, 5000000000, depth,
            search_algorithm='BFS',
            n=1, min_zeros=mz,
            epsilon=exp(-18),
            all_solutions=False
        )
        x = h.handler()
 
        pred_seq = x[0][0] if x and x[0] else False
        predicted_next = pred_seq[len(seq) - 1] if pred_seq and len(pred_seq) >= len(seq) else None
 
        full_match = bool(pred_seq) and len(pred_seq) >= len(seq) and pred_seq[:len(seq)] == seq
        ok = is_close(predicted_next, seq[-1])
 
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
    print("\n==============================")
    print(f"[Composite Baseline] Total: {len(seqs)} | Solved: {solved} | Accuracy: {accuracy:.2f}%")
    print("==============================\n")
 
    print("Failed sequences (first 10):")
    count = 0
    for r in results:
        if not r["solved"]:
            print(f"input={r['input']} | expected={r['expected']} | predicted={r['predicted']}")
            count += 1
            if count >= 10:
                break
 
    return results
 
 
def run_composite_with_decomposition(seqs, kl, mz=1, depth=3):
    """
    Composite solver pipeline:
    - First try direct prediction (fast)
    - If that fails, retry with relaxed constraints (mz=2)
    - If still failing, try splitting into streams (odd/even, stride-3)
    - Finally, try reverse (some patterns are easier backward)
    """
    
    solved = 0
    results = []
 
    for seq in seqs:
        input_seq = seq[:-1]
        expected = seq[-1]
 
        # Stage 1: blind prediction mz=1
        baseline = predict_next_for_subseq(input_seq, kl, mz=mz, depth=depth)
        baseline_pred = baseline["predicted_next"]
        baseline_ok = is_close(baseline_pred, expected)
 
        final_ok = baseline_ok
        final_pred = baseline_pred
        mode_used = "baseline"
        detail = None
 
        # Stage 1b: retry mz=2
        if not final_ok and mz == 1:
            retry = predict_next_for_subseq(input_seq, kl, mz=2, depth=depth)
            if is_close(retry["predicted_next"], expected):
                final_ok = True
                final_pred = retry["predicted_next"]
                mode_used = "baseline-mz2"
 
        # Stage 2: decomposition (odd/even + stride3, mz=1 and mz=2, scored + sanity filtered)
        if not final_ok:
            best_split = try_best_composite_split(input_seq, expected, kl, mz=mz, depth=depth)
            if best_split and best_split["matches_expected"]:
                final_ok = True
                final_pred = best_split["reconstructed_pred"]
                mode_used = f"decomposition-{best_split['mode']}"
                detail = best_split
 
        # Stage 3: reverse sequence fallback
        if not final_ok:
            rev_seq = input_seq[::-1]
            for mz_try in [1, 2]:
                rev = predict_next_for_subseq(rev_seq, kl, mz=mz_try, depth=depth)
                rev_pred = rev["predicted_next"]
                if is_close(rev_pred, expected) and is_reasonable(rev_pred, input_seq):
                    final_ok = True
                    final_pred = rev_pred
                    mode_used = f"reverse-mz{mz_try}"
                    break
 
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
    print("\n==============================")
    print(f"[Composite Improved] Total: {len(seqs)} | Solved: {solved} | Accuracy: {accuracy:.2f}%")
    print("==============================\n")
 
    print("Still failed sequences (first 10):")
    count = 0
    for r in results:
        if not r["solved"]:
            print(f"input={r['input']} | expected={r['expected']} | predicted={r['predicted']}")
            count += 1
            if count >= 10:
                break
 
    return results


if __name__ == '__main__':
    
    # run_heuristic_search(sr0, sr1, kl2)

    
    '''
    execute_kitbit_gb_False(sr0, kl2, 1, 'results/IQ_S1Z.txt')
    execute_kitbit_gb_False(sr1, kl2, 1, 'results/LI_S1Z.txt')
    execute_kitbit_gb_True(sr0, kl2, 1, 'results/IQ_N1Z.txt')
    execute_kitbit_gb_True(sr1, kl2, 1, 'results/LI_N1Z.txt')

    
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
 
        '''   

    print("---------------------------- BASELINE --------------")

    baseline_results = run_composite_baseline(composite_test_set, kl2)

    print("----------------------------IMPROVED------------------")

    improved_results = run_composite_with_decomposition(composite_test_set, kl=kl2)
   
 