
import random

class PollardRho():
    def __init__(self, a, b, n, phi_n=None):
        if not phi_n:
            self.phi_n = self.phi(n)
        self.alpha = a
        self.beta = b
        self.n = n

    # terrible way to do this
    def phi(self, n):
        amount = 0
        for k in range(1, n + 1):
            if self.ExtendedEuclidianAlgorithm(n,k)[0] == 1:
                amount += 1
        return amount

    def multiply(self, e_0, e_1):
        return (e_0 * e_1) % self.n

    def f(self, x):
        S_val = x % 3
        if S_val is 0:
            return self.multiply(x, x)
        if S_val is 1:
            return self.multiply(self.alpha, x)
        return self.multiply(self.beta, x)

    def g(self, x, a):
        S_val = x % 3
        if S_val is 0:
            return (2*a) % (self.phi_n)
        if S_val is 1:
            return (a+1) % (self.phi_n)
        return a

    def h(self, x, b):
        S_val = x % 3
        if S_val is 0:
            return (2*b) % (self.phi_n)
        if S_val is 1:
            return b
        return (b+1) % (self.phi_n)

    def raiseByExponent(self, element, exp):
        if exp is 1:
            return element
        gArray = []
        gArray.append(element)
        r = self.multiply(element, element)
        gArray.append(r)

        for i in range(2, exp.bit_length()):
            r = self.multiply(r, r)
            gArray.append(r)

        lowestBit = 0
        for i in range(exp.bit_length()):
            if ((exp & (1 << i)) != 0):
                lowestBit = i
                break
        r = gArray[lowestBit]
        lowestBit += 1

        for i in range(lowestBit, exp.bit_length()):
            if ((exp & (1 << i)) != 0):
                r = self.multiply(r, gArray[i])
        return r

    def ExtendedEuclidianAlgorithm(self, a,b):
        prevx, x = 1, 0;  prevy, y = 0, 1
        while b:
            q, r = divmod(a,b)
            x, prevx = prevx - q*x, x
            y, prevy = prevy - q*y, y
            a, b = b, r
        return a, prevx, prevy

    def modInverse(self, element, mod=None):
        if not mod:
            mod = self.n
        return self.ExtendedEuclidianAlgorithm(element, mod)[1]

    def solve(self, a_i, b_i, a_2i, b_2i, n):
        x = (b_2i-b_i) % n
        if x == 0:
            return None
        y = (a_i-a_2i) % n
        gcd, x_inv, t = self.ExtendedEuclidianAlgorithm(x, n)
        if(gcd != 1):
            if(y % gcd == 0):
                x_inv = self.modInverse(x//gcd, n//gcd)
                return ((y//gcd)*x_inv) % (n//gcd)
            return None
        return (y*x_inv) % n

    # returns the exponent for logarithm or None if failure.
    def run(self):
        if self.beta == 1:
            return 0
        if self.alpha == self.beta:
            return 1
        a_i, b_i, x_i = 0, 0, 1
        a_2i, b_2i, x_2i = 0, 0, 1

        for _ in range(1, self.n):
            x_t = x_i
            x_i = self.f(x_t)
            a_i = self.g(x_t, a_i)
            b_i = self.h(x_t, b_i)

            x_2t = x_2i
            f_2t = self.f(x_2t)
            x_2i = self.f(f_2t)
            a_2i = self.g(f_2t, self.g(x_2t, a_2i))
            b_2i = self.h(f_2t, self.h(x_2t, b_2i))
            if x_i == x_2i:
                return self.solve(a_i, b_i, a_2i, b_2i, self.n-1)
        return None


if __name__ == '__main__':
    # In multiplicative group Z_{131} the smallest primitive elements are 2 and 6.

    print(PollardRho(2, 8 ,17161).run())

    # NUM_TESTS = 1000
    # generators = [2,6]
    # G = 131
    # for g in generators:
    #     print("=================("+str(g)+")=================")
    #     fails = 0
    #     for i in range(2, G):
    #         rand_num = random.randint(2, G-1)
    #         p = PollardRho(g, i, G)
    #         res = p.run()
    #         if res is None:
    #             res = "Failure"
    #             fails +=1
    #         else:
    #             res = str(res)
    #         if res != "Failure":
    #             print("Solving log_"+str(g)+"("+str(i)+")... Result: "+res)
    #     print("==============("+str(g)+"-"+str(round(fails*100/(G-2), 2))+"%)==============\n")
