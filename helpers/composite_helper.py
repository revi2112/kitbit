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


def predict_next_for_subseq(subseq, kl, mz=1, depth=3):
    """
    Blind prediction: run KitBit on `subseq` and return the predicted next element.
    Does NOT require knowing the expected answer — no leakage.
    """
    h = KitBit(
        subseq, kl, 5000000000, depth,
        search_algorithm='BFS',
        n=1, min_zeros=mz,
        epsilon=exp(-18),
        all_solutions=False
    )
    x = h.handler()
    pred_seq = x[0][0] if x and x[0] else False
    predicted_next = pred_seq[len(subseq)] if pred_seq and len(pred_seq) > len(subseq) else None
    return {
        "predicted_next": predicted_next,
        "pred_seq": pred_seq,
        "actions": x[0][2] if x and x[0] and len(x[0]) > 2 else [],
        "time": x[1] if x and len(x) > 1 else None
    }


def split_odd_even(seq):
    return [seq[::2], seq[1::2]]


def split_stride3(seq):
    return [seq[0::3], seq[1::3], seq[2::3]]


def next_stream_index(seq_len, num_streams):
    """Position seq_len in the interleaving belongs to stream seq_len % num_streams."""
    return seq_len % num_streams


#  Sanity filter — reject predictions that are wildly out of range

def is_reasonable(pred, seq):
    """
    Reject predictions that are clearly garbage:
    - None
    - Negative (sequences in our dataset are positive)
    - More than 3x the max value in the sequence
    """
    if pred is None:
        return False
    if pred < 0:
        return False
    seq_max = max(abs(v) for v in seq) if seq else 1
    if abs(pred) > 3 * seq_max:
        return False
    return True

def score_candidate(candidate, seq):
    """
    Score a decomposition candidate. Prioritises:
      - exact match with expected (highest weight)
      - all streams produced a prediction
      - number of streams predicted
      - penalty for predictions far from the sequence's value range
    """
    pred = candidate["reconstructed_pred"]
    if pred is None:
        return -1

    seq_max = max(abs(v) for v in seq) if seq else 1
    penalty = abs(pred) / (seq_max + 1)  # normalised distance penalty

    return (
        10 * int(candidate["matches_expected"])
        + 5  * int(candidate["all_predicted"])
        + candidate["streams_predicted"]
        - 0.01 * penalty
    )


def evaluate_split(parts, expected, kl, mz=1, depth=3):
    if any(len(part) < 2 for part in parts):
        return None

    sub_results = [predict_next_for_subseq(part, kl, mz=mz, depth=depth) for part in parts]

    original_len = sum(len(p) for p in parts)
    target_stream = next_stream_index(original_len, len(parts))
    target_pred = sub_results[target_stream]["predicted_next"]

    streams_predicted = sum(1 for r in sub_results if r["predicted_next"] is not None)
    target_matches = target_pred is not None and is_close(target_pred, expected)

    return {
        "parts": parts,
        "sub_results": sub_results,
        "target_stream": target_stream,
        "reconstructed_pred": target_pred,
        "streams_predicted": streams_predicted,
        "all_predicted": streams_predicted == len(parts),
        "matches_expected": target_matches,
    }


def reconstruct_next_from_split(parts, sub_results, mode):
    """Return prediction from the stream that owns the next position."""
    original_len = sum(len(p) for p in parts)
    target_stream = original_len % len(parts)
    preds = [r["predicted_next"] for r in sub_results]
    return preds[target_stream] if preds[target_stream] is not None else None


def try_best_composite_split(seq, expected, kl, mz=1, depth=3):
  
    strategies = [("odd_even", split_odd_even), ("stride3", split_stride3)]
    candidates = []

    for mz_try in [1, 2]:
        for mode_name, splitter in strategies:
            parts = splitter(seq)
            result = evaluate_split(parts, expected, kl, mz=mz_try, depth=depth)
            if result is None:
                continue
            result["mode"] = mode_name
            result["mz_used"] = mz_try
            candidates.append(result)

    if not candidates:
        return None

    candidates.sort(key=lambda x: score_candidate(x, seq), reverse=True)

    best = candidates[0]

    if not is_reasonable(best["reconstructed_pred"], seq):
        for c in candidates[1:]:
            if is_reasonable(c["reconstructed_pred"], seq):
                return c
        return None

    return best


def is_close(a, b, tol=1e-6):
    return a is not None and b is not None and abs(a - b) < tol
