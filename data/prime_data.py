from math import isqrt


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


def generate_first_n_primes(n):
    primes = []
    x = 2
    while len(primes) < n:
        if is_prime_number(x):
            primes.append(x)
        x += 1
    return primes


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


def build_prime_families():
    primes = generate_first_n_primes(250)

    families = {
        "primes": primes,
        "prime_plus_one": [x + 1 for x in primes],
        "prime_minus_one": [x - 1 for x in primes],
        "double_primes": [2 * x for x in primes],
        "prime_squares": [x * x for x in primes],
        "prime_cubes": [x ** 3 for x in primes],
        "prime_times_index": [(i + 1) * x for i, x in enumerate(primes)],
        "nth_prime_plus_n": [x + (i + 1) for i, x in enumerate(primes)],
        "nth_prime_minus_n": [x - (i + 1) for i, x in enumerate(primes)],
        "prime_gaps": prime_gaps_list(primes),
        "twin_prime_starts": twin_prime_starts(primes),
        "cousin_prime_starts": cousin_prime_starts(primes),
        "sexy_prime_starts": sexy_prime_starts(primes),
        "cumulative_prime_sum": cumulative_sum(primes),
        "cumulative_prime_product": cumulative_product(primes, limit=12),
        "alternating_prime_sign": alternating_sign(primes),
    }

    return families


def build_prime_dataset():
    families = build_prime_families()
    results = []
    seen = set()

    def add(seq, source, strategy):
        if len(seq) < 9:
            return
        key = tuple(seq[:9])
        if key not in seen:
            seen.add(key)
            results.append({
                "sequence": list(key),
                "source": source,
                "strategy": strategy
            })

    # direct families
    for name, fam in families.items():
        for start in range(0, min(40, len(fam) - 9)):
            add(fam[start:start + 9], name, "direct")

    # pair interleavings
    selected_pairs = [
        ("primes", "prime_gaps"),
        ("primes", "prime_squares"),
        ("primes", "prime_cubes"),
        ("primes", "prime_plus_one"),
        ("primes", "prime_minus_one"),
        ("primes", "double_primes"),
        ("twin_prime_starts", "prime_gaps"),
        ("cousin_prime_starts", "prime_gaps"),
        ("sexy_prime_starts", "prime_gaps"),
        ("primes", "cumulative_prime_sum"),
    ]

    for a_name, b_name in selected_pairs:
        mixed = interleave2(families[a_name], families[b_name], max_len=120)
        for start in range(0, min(35, len(mixed) - 9)):
            add(mixed[start:start + 9], f"{a_name}+{b_name}", "pair")

    # triple interleavings
    selected_triples = [
        ("primes", "prime_gaps", "prime_squares"),
        ("primes", "prime_plus_one", "prime_minus_one"),
        ("twin_prime_starts", "prime_gaps", "double_primes"),
        ("primes", "double_primes", "prime_times_index"),
        ("primes", "nth_prime_plus_n", "nth_prime_minus_n"),
    ]

    for a_name, b_name, c_name in selected_triples:
        mixed = interleave3(families[a_name], families[b_name], families[c_name], max_len=150)
        for start in range(0, min(35, len(mixed) - 9)):
            add(mixed[start:start + 9], f"{a_name}+{b_name}+{c_name}", "triple")

    return results


prime_dataset = build_prime_dataset()
prime_test_set = [item["sequence"] for item in prime_dataset]

print(f"\n[DEBUG PRIME] Dataset size: {len(prime_test_set)}")
print(f"[DEBUG PRIME] Sample: {prime_test_set[:2]}\n")