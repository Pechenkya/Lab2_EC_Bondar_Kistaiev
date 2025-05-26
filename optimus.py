import numpy as np
import math
import random
import sympy
from Crypto.Hash import SHA256
from enum import Enum

from sympy.ntheory import factorint

class SRNG_Type(Enum):
    BBS = 1
    SHA256_DRBG = 2


class BBS:
    def __init__(self, p = int('425D2B9BFDB25B9CF6C416CC6E37B59C1F', 16), 
                        q = int('D5BBB96D30086EC484EBA3D7F9CAEB07', 16), 
                        state = 0):
        self.n = p*q
        self.state = state

        if state == 0:
            self.state = random.randint(2, self.n - 1)
    
    def generate_bytes(self, n: int):
        seq = np.zeros(n, dtype=object)
        seq[0] = self.state

        for i in range(1, n):
            seq[i] = pow(seq[i - 1], 2, self.n)

        self.state = pow(seq[-1], 2, self.n)
        seq = np.array(seq % (2**8), dtype=np.uint8) 

        return seq


class SHA256_DRBG:
    def __init__(self, seed: int, state: int | None = None, cnt = 0):
        if state:
            self.v = state.to_bytes(256 // 8)
        else:
            self.v = bytes(256 // 8)
        
        self.cnt = cnt
        self.C = SHA256.new(seed.to_bytes((seed.bit_length() + 7) // 8)).digest()
    

    @staticmethod
    def _add_bytes_modulo(*b, len=32):
        t = 0
        for bs in b:
            if isinstance(bs, int):
                t = (t + bs) % pow(2, len*8)
            elif isinstance(bs, bytes):
                t = (t + (int.from_bytes(bs))) % pow(2, len*8)
            else: 
                raise 'wtf'

        return t.to_bytes(len)

    def generate_bytes(self, n: int, addin: int | None = None):
        if addin:
            w = b'\x02' + self.v + addin.to_bytes((addin.bit_length() + 7) // 8, 'little')
            v0 = SHA256_DRBG._add_bytes_modulo(self.v, w)
        else:
            v0 = self.v

        r = bytes(0)
        v = v0
        for _ in range((n // (256 // 8)) + 1):
            r += SHA256.new(v).digest()
            v = SHA256_DRBG._add_bytes_modulo(v, 1)

        H = SHA256.new(b'\x03' + v0).digest()
        self.v = SHA256_DRBG._add_bytes_modulo(v0, H, self.C, self.cnt)
        self.cnt += 1

        return np.array([x for x in r[:n]], dtype=np.uint8) 


def rand_int(bit_length: int, gen_type = SRNG_Type.SHA256_DRBG):
    if gen_type == SRNG_Type.BBS:
        gen = BBS()
    elif gen_type == SRNG_Type.SHA256_DRBG:
        gen = SHA256_DRBG(random.randint(1024, 3456787654))

    
    req_B_len = int(math.ceil(bit_length / 8))
    r = bit_length - ((req_B_len - 1) * 8)

    while True:
        seq = gen.generate_bytes(req_B_len)
        if seq[0] < pow(2, r-1):
            continue
        
        if ("{0:b}".format(seq[0])[-r]) == '0':
            seq[0] = (seq[0] + pow(2, r-1)) % pow(2, 8)

        lastB = seq[0]
        res = int(lastB) % pow(2, r)

        for b in seq[1:req_B_len]:
            res = res*pow(2, 8) + int(b)

        return res 
        



def bytes_to_num(byte_seq):
    res = 0
    for b in byte_seq:
        res = res*(2**8) + int(b)

    return res



OPTIMUS_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 
                  61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 
                  137, 139, 149, 151, 157]
R = {}

for d in OPTIMUS_PRIMES:
    R[d] = [1]
    while R[d].count(R[d][-1]) < 2:
        R[d].append((R[d][-1] * 2) % d)
    R[d].pop()

# Метод пробних ділень
def petod_drobnyx_mylen(num):
    b = bin(num)[:1:-1]
    
    if b[0] == '0':
        return 2
    
    for d in OPTIMUS_PRIMES[1::]:
        sum = 0
        for i in range(len(b)):
            sum += int(b[i]) * R[d][i % len(R[d])]
            sum %= d
        
        if sum == 0:
            return d

    return 1
    
# Ймовірнісний алгоритм Міллера-Рабіна та загальний алгоритм для знаходження простих чисел
def miller_rabin(num, base):
    i = 1
    while (num - 1) % (2 ** i) == 0:
        i += 1

    k = i - 1
    d = (num - 1) // (2 ** k)

    a_d = pow(base, d, num)

    if a_d == 1:
        return True
    
    a_d2i = a_d
    for j in range(k):
        if a_d2i == (num - 1):
            return True
        
        a_d2i = (a_d2i ** 2) % num

    return False


def check_prime(num, error_prob = 0.0001):
    if petod_drobnyx_mylen(num) != 1:
        return False

    t = int(math.ceil(math.log(1 / error_prob, 4)))
    s = 0
    for _ in range(t):
        a = random.randrange(3, num + 1)
        s += int(miller_rabin(num, a))

    return s > (t / 2)

# Генератор простих чисел
def generate_prime_MR(bit_length: int, gen_type = SRNG_Type.SHA256_DRBG, excl = []):
    if gen_type == SRNG_Type.BBS:
        gen = BBS()
    elif gen_type == SRNG_Type.SHA256_DRBG:
        gen = SHA256_DRBG(random.randint(1024, 345604))
    
    req_B_len = int(math.ceil(bit_length / 8))
    r = bit_length - ((req_B_len - 1) * 8)

    while True:
        seq = gen.generate_bytes(req_B_len)
        if seq[0] < pow(2, r-1):
            continue
        
        if ("{0:b}".format(seq[0])[-r]) == '0':
            seq[0] = (seq[0] + pow(2, r-1)) % pow(2, 8)

        lastB = seq[0]
        res = int(lastB) % pow(2, r)

        for b in seq[1:req_B_len]:
            res = res*pow(2, 8) + int(b)

        if check_prime(res) and (res not in excl):
            return res
    


# Генератор сильнопростих чисел
def generate_safe_prime(len: int, gen_type = SRNG_Type.SHA256_DRBG, excl = []):
    if gen_type == SRNG_Type.BBS:
        gen = BBS()
    elif gen_type == SRNG_Type.SHA256_DRBG:
        gen = SHA256_DRBG(random.randint(1024, 345678))

    while True:
        seq = gen.generate_bytes(len // 8)
        if seq[0] < 128:
            continue
        
        p = bytes_to_num(seq)
        if not check_prime(p) or (p in excl):
            continue

        q = (p - 1) // 2
        if check_prime(q):
            return p
        
# Генератор блум простих чисел
def generate_blum_prime(len: int, excl = []):
    while True:
        p = generate_prime_MR(len, excl)
        if p % 4 == 3:
            return p

def generate_prime_maurer(len: int):
    if len < 64:
        pr = generate_prime_MR(len)
        return pr, list(factorint(pr).keys())
    
    while True:
        F_prime_sizes, R_size = _GenPrimeFactorSizes(len)
        F_factors = []
        F = 1
        for fs in F_prime_sizes:
            fp, _ = generate_prime_maurer(fs)
            F_factors.append(fp)
            F *= fp

        for R_tries in range(pow(R_size, 2) // 15 + 5*R_size + 2500):
            R, R_factors = _GenCompositeKnownFactors(R_size)
            if math.gcd(2*R, F) != 1:
                continue

            N = 2*F*R + 1
            
            if _CheckPrimePocklington(N, F_factors):
                return N, list(set([2] + F_factors + R_factors))
            

def _GenCompositeKnownFactors(len):
    if len <= 32:
        F = rand_int(len)
        F_factors = list(factorint(F).keys())

    else:
        tmp1, tmp2 = _GenPrimeFactorSizes(len)
        prime_sizes = tmp1 + [tmp2]

        F_factors = []
        F = 1
        for fs in prime_sizes:
            fp, _ = generate_prime_maurer(fs)
            F_factors.append(fp)
            F *= fp

    return F, F_factors

def _CheckPrimePocklington(N, F_factors):
    if petod_drobnyx_mylen(N) != 1:
        return False
    
    for q in F_factors:
        a = 2
        for base_tries in range(1000):
            if pow(a, N-1, N) != 1:
                return False
            
            if math.gcd(pow(a, (N-1) // q, N) - 1, N) == 1:
                break

            a = random.randint(3, N-2)

        if base_tries == 999:
            return False
        
    return True

def _GenPrimeFactorSizes(n_size: int):
    sizes = [0]
    while any(u < 16 for u in sizes):
        sizes = []
        while n_size - sum(sizes) > 48:
            s = random.randint(16, n_size - sum(sizes))
            sizes.append(s)
            sizes.sort(reverse=True)

        sizes.append(n_size - sum(sizes))
    

    return sizes[0:-1], sizes[-1]



def CRT(a, n):
    n_prod = math.prod(n)
    N = [n_prod // n_i for n_i in n]
    M = [pow(n_prod // n_i, -1, n_i) for n_i in n]

    return sum([a[i]*M[i]*N[i] for i in range(0, len(a))]) % n_prod


def sqrt_modp(a, p):
    a = a % p

    if sympy.jacobi_symbol(a, p) != 1:
        raise RuntimeError(f"error a = {a}, p = {p}")

    if p % 4 == 3:
        # print("4k + 3")
        sq_a = pow(a, (p + 1) // 4, p)
        return [sq_a, p - sq_a]
    
    if p % 8 == 5:
        # print("8k + 5")
        k = (p - 5) // 8
        if pow(a, 2*k + 1, p) == 1:
            sq_a = pow(a, k + 1, p)
        else:
            sq_a = (pow(a, k + 1, p) * pow(2, 2*k + 1, p)) % p

        return [sq_a, p - sq_a]
    
    if p % 8 == 1:
        # print("8k + 1")
        b = 2
        while sympy.jacobi_symbol(b, p) != -1:
            b = random.randrange(3, p - 1)

        t_a = (p - 1) // 2
        t_b = 0

        while t_a % 2 == 0:
            if (pow(a, t_a, p) * pow(b, t_b, p)) % p  == p - 1:
                t_b += (p - 1) // 2

            t_a = t_a // 2
            t_b = t_b // 2

        if (pow(a, t_a, p) * pow(b, t_b, p)) % p  == p - 1:
                t_b += (p - 1) // 2

        sq_a = (pow(a, (t_a + 1) // 2, p) * pow(b, t_b // 2, p)) % p
        return [sq_a, p - sq_a]
    

def sqrt_modpq(a, p, q):
    a = a % (p*q)

    x1, x2 = sqrt_modp(a, p)
    x3, x4 = sqrt_modp(a, q)

    sqa_1 = CRT([x1, x3], [p, q])
    sqa_2 = CRT([x1, x4], [p, q])
    sqa_3 = CRT([x2, x3], [p, q])
    sqa_4 = CRT([x2, x4], [p, q])

    return [sqa_1, sqa_2, sqa_3, sqa_4]


def check_generator(g, mod, orders):
    for qi in orders:
        if pow(g, (mod-1) // qi, mod) == 1:
            return False
        
    return True