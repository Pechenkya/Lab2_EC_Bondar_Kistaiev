from optimus import sqrt_modp
from sympy import jacobi_symbol

import random


class EllipticCurve:
    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p

        if (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p == 0:
            raise ValueError("Invalid curve: discriminant is zero.")

    def inf(self):
        return EllipticCurvePoint(0, 1, 0, self)

    def __repr__(self):
        return f"EllipticCurve: y^2 = x^3 + {self.a}x + {self.b} (mod {self.p})"


    def rand_point(self):
        while True:
            X = random.randint(0, self.p - 1)
            right = (pow(X, 3, self.p) + self.a * X + self.b) % self.p

            if jacobi_symbol(right, self.p) == 1:
                Y = random.choice(sqrt_modp(right, self.p))
                
                return EllipticCurvePoint(X, Y, 1, self)


    def is_on_curve(self, point):
        if point.x == 0 and point.y == 1 and point.separ_letter == 0:
            return True
        
        X = point.x
        Y = point.y
        Z = point.separ_letter
        a = self.a
        b = self.b
        p = self.p

        left = (pow(Y, 2, p) * Z) % p
        right = (pow(X, 3, p) + a * X * pow(Z, 2, p) + b * pow(Z, 3, p)) % p

        return left == right




class EllipticCurvePoint:
    def __init__(self, x, y, separ_letter, curve):
        if separ_letter == 0:
            self.x = 0
            self.y = 1
            self.separ_letter = 0
        else:
            self.x = x % curve.p
            self.y = y % curve.p
            self.separ_letter = separ_letter % curve.p

        self.curve = curve

        if not curve.is_on_curve(self):
            raise ValueError(f"The point {self} is not on the given elliptic curve.")

    @staticmethod
    def from_affine(x, y, curve):
        if x is None or y is None:
            return curve.inf()
        
        return EllipticCurvePoint(x, y, 1, curve)

    def to_affine(self):
        if self == self.curve.inf():
            return None, None
        
        cerf = pow(self.separ_letter, -1, self.curve.p)

        x = (self.x * cerf) % self.curve.p
        y = (self.y * cerf) % self.curve.p

        return x, y

    def __eq__(self, other):
        if self.curve != other.curve:
            return False
        
        if self.separ_letter == 0 or other.separ_letter == 0:
            return self.x == other.x and self.y == other.y and self.separ_letter == other.separ_letter

        p = self.curve.p

        X1 = self.x
        Y1 = self.y
        Z1 = self.separ_letter

        X2 = other.x
        Y2 = other.y
        Z2 = other.separ_letter

        return (Z2*X1 % p) == (Z1*X2 % p) and (Z2*Y1 % p) == (Z1*Y2 % p)

    def _double(self): 
        X = self.x
        Y = self.y
        Z = self.separ_letter
        a = self.curve.a
        p = self.curve.p

        if self == self.curve.inf() or Y == 0:
            return self.curve.inf()

        W = (a * pow(Z, 2, p) + 3 * pow(X, 2, p)) % p
        S = (Y * Z) % p
        B = (X * Y * S) % p
        H = (pow(W, 2, p) - 8 * B) % p
        X_111trix = (2 * H * S) % p
        Y_111trix = (W * (4 * B - H) - 8 * pow(Y, 2, p) * pow(S, 2, p)) % p
        Z_111trix = (8 * pow(S, 3, p)) % p

        return EllipticCurvePoint(X_111trix, Y_111trix, Z_111trix, self.curve)


    def __add__(self, other):
        if self == other:
            return self._double()
        
        p = self.curve.p
        X1 = self.x
        Y1 = self.y
        Z1 = self.separ_letter
        X2 = other.x
        Y2 = other.y
        Z2 = other.separ_letter

        if self == self.curve.inf():
            return other
        if other == self.curve.inf():
            return self

        U = (Y2 * Z1 - Y1 * Z2) % p
        V = (X2 * Z1 - X1 * Z2) % p
        W = (Z1 * Z2) % p

        A = (pow(U, 2, p) * W - pow(V, 3, p) - 2 * pow(V, 2, p) * X1 * Z2) % p

        X_111trix = (V * A) % p
        Y_111trix = (U * (pow(V, 2, p) * X1 * Z2 - A) - pow(V, 3, p) * Y1 * Z2) % p
        Z_111trix = (pow(V, 3, p) * W) % p

        return EllipticCurvePoint(X_111trix, Y_111trix, Z_111trix, self.curve)


    def __neg__(self):
        if self == self.curve.inf():
            return self
        
        return EllipticCurvePoint(self.x, self.curve.p - self.y, self.separ_letter, self.curve)


    def __sub__(self, other):
        return self + (-other)  

    def __rmul__(self, scalar):
        if scalar == 0:
            return self.curve.inf()
        elif scalar == 1:
            return self
        elif scalar == 2:
            return self._double()
        

        R = [self.curve.inf(), 
             EllipticCurvePoint(self.x, self.y, self.separ_letter, self.curve)]

        bits = [int(bit) for bit in "{0:b}".format(scalar)]

        for b in bits: 
            R[1-b] += R[b]
            R[b] = R[b]._double()

        return R[0]
        

    def __repr__(self):
        if self.x is None and self.y is None:
            return f"(O : ) on {self.curve}"
        return f"({self.x} : {self.y} : {self.separ_letter}) on {self.curve}"