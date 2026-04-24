from math import exp

from helpers.generic import KitBit


def is_prime_number(n):
    if type(n) != int or n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    i = 3
    while i * i <= n:
        if n % i == 0:
            return False
        i += 2
    return True


def next_prime_number(n):
    x = int(n) + 1
    while not is_prime_number(x):
        x += 1
    return x


def generate_first_n_primes(n):
    primes = []
    x = 2
    while len(primes) < n:
        if is_prime_number(x):
            primes.append(x)
        x += 1
    return primes


PRIME_CACHE = generate_first_n_primes(400)


def prime_gaps_list(primes):
    return [primes[i + 1] - primes[i] for i in range(len(primes) - 1)]


def twin_prime_starts(primes):
    return [primes[i] for i in range(len(primes) - 1) if primes[i + 1] - primes[i] == 2]


def cousin_prime_starts(primes):
    return [primes[i] for i in range(len(primes) - 1) if primes[i + 1] - primes[i] == 4]


def sexy_prime_starts(primes):
    return [primes[i] for i in range(len(primes) - 1) if primes[i + 1] - primes[i] == 6]


def cumulative_sum(seq):
    out = []
    s = 0
    for x in seq:
        s += x
        out.append(s)
    return out


def cumulative_product(seq, limit=12):
    out = []
    p = 1
    for x in seq[:limit]:
        p *= x
        out.append(p)
    return out


def alternating_sign(seq):
    return [x if i % 2 == 0 else -x for i, x in enumerate(seq)]


def build_prime_families():
    p = PRIME_CACHE
    fam = {}

    fam["primes"] = p
    fam["prime_plus_one"] = [x + 1 for x in p]
    fam["prime_minus_one"] = [x - 1 for x in p]
    fam["double_primes"] = [2 * x for x in p]
    fam["prime_squares"] = [x * x for x in p]
    fam["prime_cubes"] = [x ** 3 for x in p]
    fam["prime_times_index"] = [(i + 1) * x for i, x in enumerate(p)]
    fam["nth_prime_plus_n"] = [x + (i + 1) for i, x in enumerate(p)]
    fam["nth_prime_minus_n"] = [x - (i + 1) for i, x in enumerate(p)]
    fam["prime_gaps"] = prime_gaps_list(p)
    fam["twin_prime_starts"] = twin_prime_starts(p)
    fam["cousin_prime_starts"] = cousin_prime_starts(p)
    fam["sexy_prime_starts"] = sexy_prime_starts(p)
    fam["cumulative_prime_sum"] = cumulative_sum(p)
    fam["cumulative_prime_product"] = cumulative_product(p, limit=12)
    fam["alternating_prime_sign"] = alternating_sign(p)

    return fam


PRIME_FAMILIES = build_prime_families()


def solve_subsequence_with_kitbit(subseq, kl, mz=1, depth=2):
    h = KitBit(
        subseq[:-1], kl, 100000, depth,
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


def predict_next_for_subseq(subseq, kl, mz=1, depth=2):
    h = KitBit(
        subseq, kl, 100000, depth,
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


def first_differences(seq):
    return [seq[i + 1] - seq[i] for i in range(len(seq) - 1)]


def second_differences(seq):
    fd = first_differences(seq)
    return [fd[i + 1] - fd[i] for i in range(len(fd) - 1)]


def is_prime_like_sequence(seq):
    return len(seq) >= 3 and all(is_prime_number(x) for x in seq)


def match_shifted_family(seq, family):
    n = len(seq)
    if n < 3 or len(family) <= n:
        return None

    for start in range(0, len(family) - n):
        if family[start:start + n] == seq:
            return family[start + n]
    return None


def interleave2(a, b, max_len):
    out = []
    n = min(len(a), len(b))
    for i in range(n):
        out.append(a[i])
        if len(out) >= max_len:
            return out[:max_len]
        out.append(b[i])
        if len(out) >= max_len:
            return out[:max_len]
    return out[:max_len]


def interleave3(a, b, c, max_len):
    out = []
    n = min(len(a), len(b), len(c))
    for i in range(n):
        out.extend([a[i], b[i], c[i]])
        if len(out) >= max_len:
            return out[:max_len]
    return out[:max_len]


def predict_from_direct_family(seq):
    for name, family in PRIME_FAMILIES.items():
        pred = match_shifted_family(seq, family)
        if pred is not None:
            return pred, f"family-{name}"
    return None, None


def predict_from_pair_interleave_family(seq):
    family_items = list(PRIME_FAMILIES.items())

    for i in range(len(family_items)):
        for j in range(len(family_items)):
            name1, fam1 = family_items[i]
            name2, fam2 = family_items[j]

            mixed = interleave2(fam1, fam2, max_len=120)
            pred = match_shifted_family(seq, mixed)
            if pred is not None:
                return pred, f"pair-{name1}-{name2}"

    return None, None


def predict_from_triple_interleave_family(seq):
    selected = [
        "primes",
        "prime_gaps",
        "prime_squares",
        "prime_plus_one",
        "prime_minus_one",
        "double_primes",
        "twin_prime_starts",
        "cumulative_prime_sum"
    ]

    chosen = [(k, PRIME_FAMILIES[k]) for k in selected if k in PRIME_FAMILIES]

    for i in range(len(chosen)):
        for j in range(len(chosen)):
            for k in range(len(chosen)):
                name1, fam1 = chosen[i]
                name2, fam2 = chosen[j]
                name3, fam3 = chosen[k]

                mixed = interleave3(fam1, fam2, fam3, max_len=150)
                pred = match_shifted_family(seq, mixed)
                if pred is not None:
                    return pred, f"triple-{name1}-{name2}-{name3}"

    return None, None


def predict_prime_direct(seq):
    if is_prime_like_sequence(seq):
        return next_prime_number(seq[-1]), "prime-direct"
    return None, None


def predict_prime_first_diff(seq):
    diffs = first_differences(seq)
    if is_prime_like_sequence(diffs):
        next_diff = next_prime_number(diffs[-1])
        return seq[-1] + next_diff, "prime-first-diff"
    return None, None


def predict_prime_second_diff(seq):
    fd = first_differences(seq)
    sd = second_differences(seq)

    if is_prime_like_sequence(sd):
        next_sd = next_prime_number(sd[-1])
        next_fd = fd[-1] + next_sd
        return seq[-1] + next_fd, "prime-second-diff"

    return None, None


def predict_prime_odd_even(seq):
    odd_part = seq[::2]
    even_part = seq[1::2]

    for name, family in PRIME_FAMILIES.items():
        odd_pred = match_shifted_family(odd_part, family)
        even_pred = match_shifted_family(even_part, family)

        if odd_pred is not None:
            if len(odd_part) > len(even_part):
                return None, None
            return odd_pred, f"odd-even-{name}-odd"

        if even_pred is not None:
            if len(even_part) > len(odd_part):
                return None, None
            return even_pred, f"odd-even-{name}-even"

    return None, None


def predict_prime_stride3(seq):
    parts = [seq[0::3], seq[1::3], seq[2::3]]
    lengths = [len(p) for p in parts]
    min_len = min(lengths)

    for name, family in PRIME_FAMILIES.items():
        for idx, part in enumerate(parts):
            pred = match_shifted_family(part, family)
            if pred is not None:
                next_branch = None
                for j, p in enumerate(parts):
                    if len(p) == min_len:
                        next_branch = j
                        break

                if next_branch == idx:
                    return pred, f"stride3-{name}-part{idx}"

    return None, None


def deep_prime_fallback(seq):
    strategies = [
        predict_from_direct_family,
        predict_from_pair_interleave_family,
        predict_from_triple_interleave_family,
        predict_prime_direct,
        predict_prime_odd_even,
        predict_prime_stride3,
        predict_prime_first_diff,
        predict_prime_second_diff,
    ]

    for fn in strategies:
        pred, mode = fn(seq)
        if pred is not None:
            return pred, mode

    return None, None


def is_close(a, b, tol=1e-6):
    return a is not None and b is not None and abs(a - b) < tol