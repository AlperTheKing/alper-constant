"""Exact algebra and certified ball-arithmetic evaluation of the constants A(h).

Parity rock-paper-scissors game G_n(h): equal choices draw, an odd difference is won by the
smaller number, a nonzero even difference by the larger one, and a win with the number m pays
2m - h.  This program treats h = 0 (the winner collects twice his number; its constant is the
Alper constant) and h = 1.  It computes the block-sum expansion coefficients exactly, certifies
80 decimals of A(h) and, for h = 0, 30 decimals of S_0(10^m) for m = 1..20, and checks the game
identities exactly.
"""

from pathlib import Path
from math import comb
import argparse
import json
import time
import platform
import flint

from flint import arb, fmpq, fmpq_poly, fmpq_mat, ctx


def threshold(k, h):
    """First n of the block k: floor((4k^3+5k)/6)+2 for h = 1, ceil((4k^3+5k)/6)+1 for h = 0."""
    if h == 1:
        return (4 * k**3 + 5 * k) // 6 + 2
    return (4 * k**3 + 5 * k + 5) // 6 + 1


def bernoulli_numbers(n):
    b = [fmpq(1)]
    for m in range(1, n + 1):
        b.append(-sum((comb(m + 1, j) * b[j] for j in range(m)), fmpq(0)) / (m + 1))
    return b


def moment_polynomials(order):
    """E[X^j], X ~ BetaBinomial(k-1, 1/2, 3/2), exactly in Q[k]."""
    k = fmpq_poly([0, 1])
    falling = [fmpq_poly([1])]
    ratios = [fmpq(1)]
    for j in range(1, order + 1):
        falling.append(falling[-1] * (k - j))
        ratios.append(ratios[-1] * fmpq(2 * j - 1, 2 * (j + 1)))
    stirling = [1]
    moments = []
    for j in range(order + 1):
        moments.append(sum((stirling[i] * ratios[i] * falling[i] for i in range(j + 1)), fmpq_poly([])))
        stirling = [0] + [(stirling[i - 1] if i - 1 < len(stirling) else 0)
                          + (i * stirling[i] if i < len(stirling) else 0)
                          for i in range(1, j + 2)]
    return moments


def block_length_and_offset(k, sign, h):
    """ell_k = t_{k+1} - t_k - 1 and a = 2 t_k - h - 2 - (4/3)k^3 for the parity sign = (-1)^k."""
    if h == 1:
        return 2 * k**2 + 2 * k + fmpq(1 - sign, 2), fmpq(5 * k, 3) + fmpq(1 + sign, 2)
    return 2 * k**2 + 2 * k + fmpq(1 + sign, 2), fmpq(5 * k, 3) + fmpq(1 - sign, 2)


def expansion(order, sign, h, moments=None):
    """F_M(z), P_J(z), and exact rational remainder numerator G(z)."""
    mmax, degree = order, order - 1
    k = fmpq_poly([0, 1])
    if h == 1:
        length = 2 * k**2 + 2 * k + fmpq(1 - sign, 2)
        offset = fmpq(5, 3) * k + fmpq(1 + sign, 2)
    else:
        length = 2 * k**2 + 2 * k + fmpq(1 + sign, 2)
        offset = fmpq(5, 3) * k + fmpq(1 - sign, 2)
    moments = moments if moments is not None else moment_polynomials(mmax)
    bernoulli = bernoulli_numbers(mmax)
    length_powers = [length**j for j in range(mmax + 2)]
    offset_powers = [offset**j for j in range(mmax + 1)]
    powersums = []
    centered = []
    for j in range(mmax + 1):
        powersums.append(sum((comb(j + 1, l) * bernoulli[l] * length_powers[j + 1 - l]
                              for l in range(j + 1)), fmpq_poly([])) / (j + 1))
        centered.append(sum((comb(j, l) * (-4)**l * moments[l] * offset_powers[j - l]
                             for l in range(j + 1)), fmpq_poly([])))
    f = [fmpq(0) for _ in range(3 * mmax + 4)]
    for d in range(length.degree() + 1):
        f[3 - d] -= length[d] / 2
    for m in range(1, mmax + 1):
        t = sum((comb(m, j) * 2**j * powersums[j] * centered[m - j]
                 for j in range(m + 1)), fmpq_poly([]))
        factor = -fmpq(-3, 4)**m
        for d in range(t.degree() + 1):
            p = 3 * m + 1 - d
            assert p >= 0
            f[p] += factor * t[d]
            f[p + 2] += factor * t[d] / 2
    coefficients = []
    for j in range(degree + 1):
        coefficients.append((f[j] - (coefficients[-1] if j else 0)) / 2)
    g = f.copy()
    for j, coefficient in enumerate(coefficients):
        g[j] -= 2 * coefficient
        g[j + 1] -= coefficient
    assert all(v == 0 for v in g[:degree + 1])
    leading = fmpq(-3 * sign, 8) if h == 1 else fmpq(3 * sign, 8)
    assert coefficients[:3] == [fmpq(3, 2), fmpq(-3, 4), leading]
    return coefficients, g


def power_tail_bound(s, cutoff):
    return fmpq(1, cutoff**s) + fmpq(1, (s - 1) * cutoff**(s - 1))


def truncation_bound(order, cutoff, remainders):
    assert cutoff >= 8 and order >= 3
    geometric = fmpq(3, 2) * 4**(order + 1) / (1 - fmpq(4, cutoff)) * power_tail_bound(order, cutoff)
    polynomial = sum((max(abs(remainders[0][p]), abs(remainders[1][p])) * power_tail_bound(p, cutoff) / 2
                      for p in range(order, len(remainders[0]))), fmpq(0))
    return geometric, polynomial


def central_coefficients(n):
    c = [fmpq(1)]
    for j in range(1, n):
        c.append(c[-1] * fmpq(2 * j - 1, 2 * j))
    return c


def block_ball(k, c, h, first=None, last=None):
    """B_k (or a partial block) from the digamma formula; harmonic_j = sum_n 1/(2n - h - 4j + 2)."""
    first = threshold(k, h) if first is None else first
    last = threshold(k + 1, h) - 2 if last is None else last
    harmonic = (arb(fmpq(2 * last - h, 2)).digamma()
                - arb(fmpq(2 * first - 2 - h, 2)).digamma()) / 2
    weighted = arb(0)
    for x in range(k):
        weight = 2 * c[x] * (2 * k - 2 * x - 1) * c[k - 1 - x]
        weighted += arb(weight) * harmonic
        if x + 1 < k:
            harmonic += (arb(fmpq(1, 2 * first - h - 4 * x - 6))
                         + arb(fmpq(1, 2 * first - h - 4 * x - 4))
                         - arb(fmpq(1, 2 * last - h - 4 * x - 4))
                         - arb(fmpq(1, 2 * last - h - 4 * x - 2)))
    return (last - first + 1 - arb(fmpq(2 * k*k + 1, 3)) * weighted) / (2 * k + 1)


def exact_probability(k, n, h):
    """P_k at the largest payment Y = 2n - h: R_k = prod (Y - 4j)/(Y - 4j + 2)."""
    r = fmpq(1)
    for j in range(1, k + 1):
        r *= fmpq(2 * n - h - 4 * j, 2 * n - h - 4 * j + 2)
    return ((2 * k*k + 1) * r - 2 * (k*k - 1)) / (3 * (2*k + 1))


def verify_small(h):
    c = central_coefficients(15)
    for k in range(1, 15):
        exact = sum((exact_probability(k, n, h) for n in range(threshold(k, h), threshold(k + 1, h) - 1)), fmpq(0))
        assert block_ball(k, c, h).overlaps(arb(exact))
        if k == 1:
            assert exact == (fmpq(11239, 10395) if h == 1 else fmpq(163, 180))
        for n in [threshold(k, h), threshold(k + 1, h) - 2]:
            r = (3 * (2*k+1)*exact_probability(k, n, h) + 2*(k*k-1))/(2*k*k+1)
            partial = 2 * sum(((2*q+1)*c[q]*c[k-1-q] / (2*n-h-4*k+4*q+2)
                               for q in range(k)), fmpq(0))
            assert 1-r == partial
    print('Exact product, partial fractions and digamma checks passed.', flush=True)


def payoff(i, j, h):
    if i == j:
        return 0
    winner = min(i, j) if (i-j) % 2 else max(i, j)
    return (2*winner-h) * (1 if winner == i else -1)


def verify_game(h, limit=200):
    for n in range(3, limit + 1):
        k = 1
        while threshold(k + 1, h) <= n:
            k += 1
        last = n-1 if n == threshold(k+1, h)-1 else n
        first, size = last-2*k, 2*k+1
        matrix = fmpq_mat([[payoff(first+i, first+j, h) for j in range(size)]
                          for i in range(size-1)] + [[1]*size])
        vector = matrix.solve(fmpq_mat([[0]]*(size-1) + [[1]]))
        assert all(vector[j, 0] > 0 for j in range(size))
        for i in range(1, n+1):
            f = sum((payoff(i, first+j, h)*vector[j, 0] for j in range(size)), fmpq(0))
            assert f == 0 if first <= i <= last else f < 0
        assert vector[size-1, 0] == exact_probability(k, last, h)
        r = (3*(2*k+1)*vector[size-1, 0] + 2*(k*k-1))/(2*k*k+1)
        lower = sum((payoff(first-1, first+j, h)*vector[j, 0] for j in range(size)), fmpq(0))
        upper = sum((payoff(last+1, first+j, h)*vector[j, 0] for j in range(size)), fmpq(0))
        assert lower == fmpq(2*k+3, 2*k+1)*(2*last-h)*r*exact_probability(k+1, last+1, h)
        assert upper == -fmpq(2*k+3, 2*k+1)*(2*last-h+2)*exact_probability(k+1, last+2, h)
    print('Exact full-game equilibria and both neighbor identities verified through n =', limit, flush=True)


def verify_moments(order=20, max_k=35):
    moments = moment_polynomials(order)
    c = central_coefficients(max_k)
    for k in range(1, max_k+1):
        sums = [fmpq(0)]*(order+1)
        for x in range(k):
            weight = c[x]*(2*k-2*x-1)*c[k-1-x]/k
            power = 1
            for j in range(order+1):
                sums[j] += weight*power
                power *= x
        assert all(moments[j](k) == sums[j] for j in range(order+1))
    print('Factorial-moment polynomials independently checked against all weights.', flush=True)


def tail_ball(cutoff, coeff, remain, order):
    tail = arb(0)
    for parity in [0, 1]:
        first = cutoff + (parity-cutoff) % 2
        for j in range(2, order):
            tail += arb(coeff[parity][j])*arb(j).zeta(arb(fmpq(first, 2)))/2**j
    geometric, polynomial = truncation_bound(order, cutoff, remain)
    return tail + arb(0, arb(geometric+polynomial).upper())


def partial_block_ball(k, length, order, moments, h):
    if length == 0:
        return arb(0)
    sign = (-1)**k
    _, offset = block_length_and_offset(k, sign, h)
    bernoulli = bernoulli_numbers(order)
    sums = [sum((comb(j+1, l)*bernoulli[l]*length**(j+1-l)
                 for l in range(j+1)), fmpq(0))/(j+1) for j in range(order+1)]
    shifted = [sum((comb(j, l)*(-4)**l*moments[l](k)*offset**(j-l)
                    for l in range(j+1)), fmpq(0)) for j in range(order+1)]
    numerator = arb(-fmpq(length, 2*k*k))
    for m in range(1, order+1):
        tm = sum((comb(m, j)*2**j*sums[j]*shifted[m-j] for j in range(m+1)), fmpq(0))
        numerator -= arb((1+fmpq(1, 2*k*k))*fmpq(-3, 4*k**3)**m*tm)
    result = numerator/(2*k+1)
    error = fmpq(3*k, 2)*fmpq(4, k)**(order+1)/(1-fmpq(4, k))
    return result+arb(0, arb(error).upper())


def partial_block_centered_ball(k, length, h, order=12):
    """Independent pole-centered expansion with exact Hurwitz-zeta sums."""
    first = threshold(k, h)
    a = 2*first-h-k-1       # mean of 2n - h - 2 - 4X at n = first
    b = a+2*length
    spread = k-1
    moments = moment_polynomials(order)
    weighted = arb(0)
    for j in range(order+1):
        moment = sum((comb(j, l)*4**l*(-spread)**(j-l)*moments[l](k)
                      for l in range(j+1)), fmpq(0))
        if j == 0:
            harmonic = (arb(fmpq(b, 2)).digamma()-arb(fmpq(a, 2)).digamma())/2
        else:
            harmonic = (arb(j+1).zeta(arb(fmpq(a, 2)))
                        -arb(j+1).zeta(arb(fmpq(b, 2))))/2**(j+1)
        weighted += arb(moment)*harmonic
    coefficient = fmpq(2*k*(2*k*k+1), 3*(2*k+1))
    value = arb(fmpq(length, 2*k+1))-arb(coefficient)*weighted
    rho = fmpq(3*spread, a)
    error = coefficient*fmpq(length, a)*rho**(order+1)/(1-rho)
    return value+arb(0, arb(error).upper())


def rounded_decimal(ball, places=80):
    scaled = (ball*10**places+arb(fmpq(1, 2))).floor().unique_fmpz()
    if scaled is None:
        raise ArithmeticError('The enclosure does not prove the requested decimal rounding.')
    digits = str(int(scaled))
    digits = digits.zfill(places+1)
    return digits[:-places]+'.'+digits[-places:]


def sum_at(n, total, coeff, remain, order, moments, h):
    low, high = 1, 2
    while threshold(high, h) <= n:
        high *= 2
    while high-low > 1:
        mid = (low+high)//2
        if threshold(mid, h) <= n:
            low = mid
        else:
            high = mid
    k = low
    assert k >= 8
    length = min(n, threshold(k+1, h)-2)-threshold(k, h)+1
    harmonic = arb(k).digamma()-arb(1).digamma()
    partial = partial_block_ball(k, length, order, moments, h)
    other = partial_block_centered_ball(k, length, h)
    assert partial.overlaps(other)
    result = (arb(fmpq(3*(k-1), 2))-3*harmonic/4+total
              -tail_ball(k, coeff, remain, order)
              +partial)
    return k, result


def compute(h, cutoff, order, digits, output, powers_table=True):
    if cutoff < 8 or order < 3 or digits < 50:
        raise ValueError('Require cutoff >= 8, order >= 3 and digits >= 50.')
    ctx.dps = digits
    start = time.perf_counter()
    verify_small(h)
    verify_game(h)
    verify_moments()
    moments = moment_polynomials(order)
    coeff, remain = [], []
    for sign in [1, -1]:
        a, g = expansion(order, sign, h, moments)
        coeff.append(a)
        remain.append(g)
    check_c = central_coefficients(cutoff+2)
    for k in [cutoff, cutoff+1]:
        parity = k % 2
        approximate = sum((arb(coeff[parity][j])/k**j for j in range(order)), arb(0))
        bound = fmpq(3*k, 2)*fmpq(4, k)**(order+1)/(1-fmpq(4, k))
        bound += sum((abs(remain[parity][p])*fmpq(1, k**p)/2
                      for p in range(order, len(remain[parity]))), fmpq(0))
        assert (approximate+arb(0, arb(bound).upper())).overlaps(block_ball(k, check_c, h))
    print('Exact polynomial expansions constructed.', flush=True)
    for j in range(min(11, order)):
        print('coefficient', j, 'even', coeff[0][j], 'odd', coeff[1][j], flush=True)
    geometric, polynomial = truncation_bound(order, cutoff, remain)
    error = arb(geometric + polynomial)
    c = central_coefficients(cutoff)
    finite = sum((block_ball(k, c, h) - arb(fmpq(3, 2)) + arb(fmpq(3, 4*k))
                  for k in range(1, cutoff)), arb(0))
    tail = arb(0)
    for parity in [0, 1]:
        first = cutoff + (parity - cutoff) % 2
        for j in range(2, order):
            tail += arb(coeff[parity][j]) * arb(j).zeta(arb(fmpq(first, 2))) / 2**j
    total = finite + tail + arb(0, error.upper())
    cprime = arb(fmpq(3, 2)) - 3 * arb(1).digamma() / 4 - total
    constant = cprime + arb(fmpq(3, 2)).log() / 4
    # leading part: e_k ~ sigma (-1)^(k+1) 3/(8k^2), sigma = +1 for h = 1 and -1 for h = 0
    sigma = 1 if h == 1 else -1
    irregular = total - sigma * arb.pi()**2 / 32
    powers = {}
    for m in (range(1, 21) if powers_table else []):
        n = 10**m
        if m <= 9:  # every block exactly, by the digamma formula
            kk = 1
            while threshold(kk + 1, h) <= n:
                kk += 1
            cc = central_coefficients(kk + 2)
            ball = sum((block_ball(j, cc, h) for j in range(1, kk)), arb(0))
            ball += block_ball(kk, cc, h, threshold(kk, h), min(n, threshold(kk + 1, h) - 2))
        else:
            ball = sum_at(n, total, coeff, remain, order, moments, h)[1]
        powers[m] = ball
    result = {
        'h': h,
        'python_version': platform.python_version(), 'python_flint_version': flint.__version__,
        'cutoff': cutoff, 'order': order, 'digits': digits,
        'geometric_bound': str(arb(geometric)), 'polynomial_bound': str(arb(polynomial)),
        'finite_sum': str(finite), 'tail': str(tail),
        'sum_e': str(total), 'C_prime': str(cprime), 'A': str(constant),
        'remainder_after_pi2_32': str(irregular), 'pi2_32_sign': sigma,
        'radius_sum_e': str(total.rad()),
        'coefficients_even': [str(v) for v in coeff[0]],
        'coefficients_odd': [str(v) for v in coeff[1]],
        'seconds': time.perf_counter() - start,
    }
    if powers_table:
        result['block_at_1e20'], sum_1e20 = sum_at(10**20, total, coeff, remain, order, moments, h)
        result['S_1e20'] = str(sum_1e20)
        result['S_powers_of_ten'] = {str(m): str(b) for m, b in powers.items()}
        result['partial_block_independent_ball_check'] = True
    if platform.system() == 'Windows':
        import ctypes
        count = ctypes.windll.kernel32.GetActiveProcessorCount
        count.argtypes = [ctypes.c_ushort]
        count.restype = ctypes.c_uint
        result['logical_processors_all_groups'] = count(0xffff)
    result['threads_used'] = 1
    if total.rad() < arb('1e-85'):
        result['certified_80_decimals'] = {name: rounded_decimal(value)
            for name, value in [('sum_e', total), ('C_prime', cprime), ('A', constant)]}
        if powers_table:
            result['certified_80_decimals']['S_1e20'] = rounded_decimal(sum_1e20)
            result['certified_30_decimals_S_10m'] = {str(m): rounded_decimal(b, 30) for m, b in powers.items()}
    Path(output).write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    for key in ['geometric_bound', 'polynomial_bound', 'sum_e', 'C_prime', 'A', 'remainder_after_pi2_32', 'radius_sum_e', 'S_1e20', 'seconds']:
        if key in result:
            print(key + ':', result[key], flush=True)
    return result


def research_checks(certificate_path, output):
    import mpmath as mp
    import sympy as sp
    data = json.loads(Path(certificate_path).read_text(encoding='utf-8'))
    h = data['h']
    mp.mp.dps = 130
    exact_value = lambda q: mp.mpf(str(q.numerator))/int(q.denominator)
    cutoff = 80
    c = central_coefficients(cutoff)
    independent = mp.mpf(0)
    half = mp.mpf(h)/2
    # Deliberately use a separate digamma pair for every pole (no moving window).
    for k in range(1, cutoff):
        a, b = threshold(k, h), threshold(k+1, h)-2
        weighted = mp.mpf(0)
        for x in range(k):
            weight = 2*c[x]*(2*k-2*x-1)*c[k-1-x]
            hx = (mp.digamma(mp.mpf(b)-half-2*x)-mp.digamma(mp.mpf(a)-half-1-2*x))/2
            weighted += exact_value(weight)*hx
        independent += ((b-a+1)-mp.mpf(2*k*k+1)*weighted/3)/(2*k+1)-mp.mpf('1.5')+mp.mpf(3)/(4*k)
    for parity, key in enumerate(['coefficients_even', 'coefficients_odd']):
        a = mp.mpf(cutoff+(parity-cutoff) % 2)/2
        for j in range(2, len(data[key])):
            independent += exact_value(fmpq(data[key][j]))*mp.zeta(j, a)/2**j
    target = mp.mpf(data['sum_e'].split()[0].lstrip('['))
    difference = abs(independent-target)
    print('Independent mpmath sum_e:', mp.nstr(independent, 90), flush=True)
    print('Difference from certified midpoint:', mp.nstr(difference, 8), flush=True)
    assert difference < mp.mpf('1e-70')

    k, n, r, hs = sp.symbols('k n r h')
    length = 2*k+1
    y = n-hs/2
    p = ((2*k*k+1)*r-2*(k*k-1))/(3*length)
    q = (2*y*(r-1)+2*k)/3+1/length
    nextp = lambda rr: ((2*(k+1)**2+1)*rr-2*((k+1)**2-1))/(3*(2*k+3))
    lower = 2*y*p-length*(1+q)
    upper = length+1/length-2*q-2*y*p
    assert sp.cancel(lower-(2*k+3)/length*2*y*r*nextp((y-2*k-1)/(y*r))) == 0
    assert sp.cancel(upper+(2*k+3)/length*(2*y+2)*nextp(y*r/(y+1))) == 0
    gosper = {}
    for kk in [1, 2, 3, 4]:
        rr = sp.prod((2*n-h-4*j)/(2*n-h-4*j+2) for j in range(1, kk+1))
        pp = sp.cancel(((2*kk*kk+1)*rr-2*(kk*kk-1))/(3*(2*kk+1)))
        gosper[str(kk)] = str(sp.concrete.gosper.gosper_sum(pp, n))

    mp.mp.dps = 120
    values = {
        'sum_e': target,
        'C_prime': mp.mpf('1.5')+mp.mpf('0.75')*mp.euler-target,
        'A': mp.mpf('1.5')+mp.mpf('0.75')*mp.euler-target+mp.log(mp.mpf('1.5'))/4,
    }
    bases = {
        'basic': [mp.mpf(1), mp.euler, mp.log(2), mp.log(3), mp.pi**2, mp.zeta(3)],
        'extended': [mp.mpf(1), mp.euler, mp.log(2), mp.log(3), mp.pi, mp.pi**2,
                     mp.pi**4, mp.zeta(3), mp.zeta(5), mp.catalan, mp.log(mp.pi),
                     mp.loggamma(mp.mpf(1)/4), mp.loggamma(mp.mpf(1)/3), mp.sqrt(2),
                     mp.pi*mp.sqrt(3), mp.log(2)**2, mp.log(3)**2,
                     mp.log(2)*mp.log(3), mp.pi*mp.log(2)]
    }
    searches = {}
    for name, value in values.items():
        for basis, numbers in bases.items():
            relation = mp.pslq(mp.matrix([value]+numbers), tol=mp.mpf('1e-85'), maxcoeff=1000, maxsteps=10000)
            searches[name+'_'+basis] = relation
            print('PSLQ', name, basis, relation, flush=True)
            if relation:
                raise RuntimeError('A PSLQ candidate needs a fresh higher-precision validation.')
    result = {'h': h, 'independent_sum_e': mp.nstr(independent, 100),
              'difference_from_certificate': mp.nstr(difference, 15),
              'neighbor_identities_symbolic': True, 'gosper_in_n': gosper,
              'pslq': searches, 'pslq_dps': 120, 'pslq_tolerance': '1e-85',
              'pslq_maxcoeff': 1000, 'pslq_maxsteps': 10000}
    Path(output).write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--h', type=int, choices=[0, 1], default=0, help='payment 2m - h; h = 0 gives the Alper constant')
    parser.add_argument('--cutoff', type=int, default=192)
    parser.add_argument('--order', type=int, default=60)
    parser.add_argument('--digits', type=int, default=160)
    parser.add_argument('--output')
    parser.add_argument('--powers', choices=['auto', 'yes', 'no'], default='auto',
                        help='S_h(10^m) table, m <= 20; auto = only for h = 0 (the values for h = 1 are a Project Euler answer)')
    parser.add_argument('--research-checks', metavar='CERTIFICATE')
    args = parser.parse_args()
    if args.research_checks:
        research_checks(args.research_checks, args.output or 'research_checks.json')
    else:
        powers_table = args.powers == 'yes' or (args.powers == 'auto' and args.h == 0)
        compute(args.h, args.cutoff, args.order, args.digits, args.output or 'certificate_h%d.json' % args.h, powers_table)
