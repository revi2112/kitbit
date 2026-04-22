from math import exp

from helpers.generic import KitBit


def solve_subsequence_with_kitbit(subseq, kl, mz=1, depth=3):
    h = KitBit(
        subseq[:-1], kl, 5000000000, depth,
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