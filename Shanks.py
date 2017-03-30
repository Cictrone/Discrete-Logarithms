import math

# written for Z_{n^2}[x]
class Shanks():
    def __init__(self, a, b, n):
        self.a = a
        self.b = b
        self.n = n
    # written for Z_{n^2}[x] modulo x^2 +1, a and b are tuples modulo n
    def multiply(self, a, b):
        h = (a[0]*b[1] + a[1]*b[0]) % self.n
        k = (a[1]*b[1] - a[0]*b[0]) % self.n
        return (h,k)

    # Will find log_alpha(b) with Shanks' Algorithm
    def run(self):
        order = self.n*self.n - 1
        m = math.ceil(math.sqrt(order))
        L_1 = []
        L_2 = []

        alpha_inv = self.raiseByExponent(self.a, order-1)

        for j in range(m):
            L_1.append((j, self.raiseByExponent(self.a, m*j)))

            alpha_neg_j = self.raiseByExponent(alpha_inv, j)
            L_2.append((j, self.multiply(alpha_neg_j, self.b)))

        L_1.sort(key=lambda x: x[1][0]*self.n + x[1][1])
        L_2.sort(key=lambda x: x[1][0]*self.n + x[1][1])
        for first in L_1:
            for second in L_2:
                if(first[1] == second[1]):
                    return (m*first[0] + second[0]) % order
        return None

    def raiseByExponent(self,  element, exp):
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


if __name__ == '__main__':
    p = 131
    alpha = (1,3) # x-3
    beta = (1, 101) # x-101
    s = Shanks(alpha, beta, p)
    ans = s.run()
    print(ans)

    print(s.raiseByExponent(alpha, ans))
