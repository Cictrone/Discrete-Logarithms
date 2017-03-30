from PollardRho import PollardRho

class PohligHellman():
    def  __init__(self, a, b, n, qList, cList):
        self.a = a
        self.b = b
        self.n = n
        self.qList = qList
        self.cList = cList

    def multiply(self, e_0, e_1):
        return (e_0 * e_1) % self.n

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

    def ExtendedEuclidianAlgorithm(self, a, b):
        prevx, x = 1, 0;  prevy, y = 0, 1
        while b:
            q, r = divmod(a,b)
            x, prevx = prevx - q*x, x
            y, prevy = prevy - q*y, y
            a, b = b, r
        #    [gcd]=a[x] + b[y]
        return a, prevx, prevy

    def modInverse(self, element, p=None):
        # assumes inverse exists, i.e. gcd(e,n)=1
        if not p:
            p = self.n
        return self.ExtendedEuclidianAlgorithm(element, p)[1]

    def naiveLog(self, alpha, beta, n):
        val = alpha
        for i in range(2, n):
            val = (val*alpha) % n
            if(val == beta):
                return i
        return None

    def algorithm(self, n, alpha, beta, q, c):
        x_i = []
        for j in range(c):
            exp = n//(q**(j+1))
            delta = self.raiseByExponent(beta, exp)
            base = self.raiseByExponent(alpha, n//q)
            x_j = PollardRho(base, delta, n+1).run()
            if not x_j:
                x_j = self.naiveLog(base, delta, n+1)
            x_i.append(x_j)
            exp = self.multiply(x_j, q**j)
            inv_mult = self.raiseByExponent(alpha, exp)
            if (beta % inv_mult == 0):
                beta = beta // inv_mult
            else:
                multiplier = self.modInverse(inv_mult)
                beta = self.multiply(beta, multiplier)
        return x_i

    # each element is a tuple of remainder and modulus
    def ChineseRemainderTheorem(self, CRTComponents):
        y = 0
        for i in range(len(CRTComponents)):
            y_i = CRTComponents[i][0]
            t_i = 1
            for j in range(len(CRTComponents)):
                if i != j:
                    t_i *= CRTComponents[j][1]
            t_inv_i = self.modInverse(t_i, CRTComponents[i][1])
            y_i *= t_i
            y_i *= t_inv_i
            y += y_i
        m = 1
        for comp in CRTComponents:
            m *= comp[1]
        return y % m

    def run(self):
        CRTComponents = []
        for i in range(len(self.qList)):
            q_i = self.qList[i]
            c_i = self.cList[i]
            x_q = self.algorithm(self.n-1, self.a, self.b, q_i, c_i)
            x_i = 0
            for j in range(len(x_q)):
                x_i += self.multiply(x_q[j], q_i**j)
            CRTComponents.append((x_i, q_i**c_i))
        return self.ChineseRemainderTheorem(CRTComponents)


if __name__ == '__main__':
    a = 2
    b = 103
    n = 131**2
    qList = [2,5,13,131]
    cList = [1,1,1,1]

    test = PohligHellman(a, b, n, qList, cList)
    print(test.run())
    # print("I got:", test.run(),"as my exponent.")
    # print("That means that", test.raiseByExponent(a, test.run()), "equals", b)
