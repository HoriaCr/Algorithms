import random
import math
import string
from collections import deque

def prime(val):
    if val == 2: 
        return 1
    if val == 1 or val % 2 == 0:
        return 0
    i = 3
    while i*i <= val:
        if val%i == 0:
            return 0
        i += 2
        
        
    return 1
    
   
#for i in range(2,30):
#   if prime(i) == 1 :
#       print i
        

firstName = 'john'
lastName = 'snowflake'

#print firstName + ' ' + lastName



def expectedNumberOfGuesses(val):
    lowBound = 1
    highBound = 100
    guessesMade = 0
    while lowBound <= highBound:
        mid = (lowBound + highBound)>>1
        guessesMade += 1
        if mid == val:
            return guessesMade
        if mid > val:
            highBound = mid - 1
        else:
            lowBound = mid + 1
            
    return -1

def game():
    name = raw_input('Hello! What is your name?\n')
    
    while True:
        guesses_made = 0
        
        number = random.randint(1, 100)
        print 'Well, {0}, I am thinking of a number between 1 and 100.'.format(name)
        
        while guesses_made < 6:
        
            guess = int(raw_input('Take a guess: '))
        
            guesses_made += 1
        
            if guess < number:
                print 'Your guess is too low.'
        
            if guess > number:
                print 'Your guess is too high.'
        
            if guess == number:
                break
        
        if guess == number:
            print 'Good job, {0}! You guessed my number in {1} guesses! /n The expected number of guesses via binary search is {2}'.format(name, guesses_made,expectedNumberOfGuesses(number))
        else:
            print 'Nope. The number I was thinking of was {0}'.format(number)
            
        inputData = raw_input( 'if you want to play again type "yes" ')
        if inputData != 'yes':
            break
        
#euler project problems in python

#the Sum of all multiples of 3 or 5 below 1000
def problem1():
    Sum = 0
    for i in range(1000):
        if i%3 == 0 or i%5 == 0:
            Sum += i
#print 'the Sum of all multiples of 3 or 5 below 1000 is {0}'.format(Sum)

#By considering the terms in the Fibonacci 
#sequence whose values do not exceed four million, find the Sum of the even-valued terms

def problem2():
    fiba = 1
    fibb = 1
    fibc = 0
    Sum = 0
    while True:
        fibc = fiba + fibb
        fiba = fibb
        fibb = fibc
        if fibc > 4000000:
            break
        
        if fibc%2 == 0:
            Sum += fibc
            
    print 'the Sum is {0}'.format(Sum)


#What is the largest prime factor of the number 600851475143 ?

def problem3(value):
    d = 2
    ret = 1
    while d*d <= value:
        if value%d == 0:
            while value%d == 0:
                value /= d
        ret = max(ret,d)
        d += 1
    ret = max(ret,value)
    
    return ret

#print 'the largest factor of {0} is : {1} '.format(600851475143,largestFactor(600851475143))

#Find the largest palindrome made from the product of two 3-digit numbers.

def palindrome(value):
    val = str(value)
    n = len(val)
    for i in range(0,n/2 + 1):
        if val[i] != val[n - i - 1]:
            return False
    return True

def problem4():
    ret = -1
    for i in range(100,1000):
        for j in range(100,1000):
            if palindrome(i*j) == True:
                ret = max(ret,i*j)
    
    return ret

#print 'the largest palindrome is : {0}'.format(largestPalindrome())

def gcd(a,b):
    return a if b == 0 else gcd(b,a%b)

def lcm(a,b):
    return a*b/gcd(a,b)


#What is the smallest positive number that is 
#evenly divisible by all of the numbers from 1 to 20?

def problem5():
    ret = 1
    for i in range(1,21):
        ret = lcm(i,ret)
    
    return ret

#print problem5()

def problem6(n):
    Sum2 = n*(n + 1)/2
    Sum2 *= Sum2
    Sum1 = n*(n + 1)*(2*n + 1)/6
    return Sum2 - Sum1

#print problem6(100)

#What is the 10 001st prime number?
def problem7(n):
    p = [0 for i in range(500000)]
    k = 0
    for i in range(2,500000):
        if p[i] == 0:
            k += 1
            if k == n:
                return i
            
            j = i + i
            while j < 500000:
                p[j] = 1
                j += i
    return 'it takes too much time'

#print 'the 10001st prime number is : {0}'.format(problem7(10001))


#Find the greatest product of five consecutive digits in the 1000-digit number.
def problem8():
    ret = 0
    num = '7316717653133062491922511967442657474235534919493496983520312774506326239578318016984801869478851843858615607891129494954595017379583319528532088055111254069874715852386305071569329096329522744304355766896648950445244523161731856403098711121722383113622298934233803081353362766142828064444866452387493035890729629049156044077239071381051585930796086670172427121883998797908792274921901699720888093776657273330010533678812202354218097512545405947522435258490771167055601360483958644670632441572215539753697817977846174064955149290862569321978468622482839722413756570560574902614079729686524145351004748216637048440319989000889524345065854122758866688116427171479924442928230863465674813919123162824586178664583591245665294765456828489128831426076900422421902267105562632111110937054421750694165896040807198403850962455444362981230987879927244284909188845801561660979191338754992005240636899125607176060588611646710940507754100225698315520005593572972571636269561882670428252483600823257530420752963450'
    zeroCode = ord('0')
    for i in range(4,len(num)):
        val = 1
        for j in range(i - 4,i + 1):
            val *= (ord(num[j]) - zeroCode)
            
        ret = max(ret,val)
    return ret

#print 'the answer is {0}'.format(problem8())

def problem9():
    for a in range(1,1000):
        b = a + 1
        while a + b < 1000:
            c = 1000 - a - b
            if a*a + b*b == c*c:
                return a*b*c
            b += 1
            
    return -1

#print 'the answer is {0}'.format(problem9())

def problem10():
    ret = 0
    p = [0 for i in range(2000000)]
    for i in range(2,2000000):
        if p[i] == 0:
            ret += i
            j = i + i
            while j < 2000000:
                p[j] = 1
                j += i
    return ret

#print 'the Sum of all the primes below two million is {0}'.format(problem10())

def getDigit(d): 
    return ord(d) - 48



def problem11():
    ret = 0
    fname = "input.txt"
    with open(fname) as f:
        lines = f.readlines()
    
    A = []
    i = 0
    for line in lines:
        j = 0
        A.append([])
        while j < 59:
            A[i].append( getDigit(line[j])*10 + getDigit(line[j + 1]))
            j += 3 
        i += 1
    dx = [-1,-1,-1,1,1,1,0,0]
    dy = [-1,1,0,-1,1,0,-1,1]
    for i in range(0,20):
        for j in range(0,20):
            for k in range(0,8):
                x = i + dx[k]
                y = j + dy[k]
                p = A[i][j]
                w = 0
                while x >= 0 and y >= 0 and x < 20 and y < 20 and w < 3:
                    p *= A[x][y]
                    x += dx[k]
                    y += dy[k]
                    w += 1
                  
                if w == 3:
                    ret = max(ret,p)
    
    return ret
            
        

#print problem11()

def getNumberOfDivisors(val):
    num = 1
    currPower = 0
    while val%2 == 0:
        val /= 2
        currPower += 1
     
    num *= (currPower + 1)
    
    d = 3
    while d*d <= val:
        currPower = 0
        while val%d == 0:
            val /= d
            currPower +=1
        
        num *= (currPower + 1)
        d += 2
            
    if val > 1:
        num *= 2
    
    return num    

def problem12():
    triangleNumber = 0
    for i in range(1000,100000):
        triangleNumber = i*(i + 1)/2
        if getNumberOfDivisors(triangleNumber) > 500:
            print i
            return triangleNumber
        
    return -1


#print problem12()

def problem13():
    fname = "input.txt"
    with open(fname) as f:
        lines = f.readlines()
  
    ret = 0
    for line in lines:
        value = 0
        for i in range(0,50):
            value = value*10 + getDigit(line[i])
        #print value
        ret += value
            
    return ret    

#print problem13()
#colla = [0 for i in range(0,10000000)]
#colla[1] = 1
def collatz(n):
    if n == 1:
        return 1
    
    return 1 + collatz( n/2 if n%2 == 0 else 3*n + 1)
    #if (colla[n] != 0):
    #   return colla[n]
    #colla[n] = 1 + collatz( n/2 if n%2 == 0 else 3*n + 1)
    #return colla[n]

def problem14():
    ret = 0
    maxLength = 0
    for i in range(1,1000000):
        chainLength = collatz(i)
        if chainLength > maxLength:
            ret = i
            maxLength = chainLength
    return ret

#print problem14()

dp = []
for i in range(0,64):
    dp.append([-1 for j in range(0,64)])
    
dp[0][0] = 1

def problem15(i,j):
    if dp[i][j] != -1:
        return dp[i][j]
    dp[i][j] = 0
    
    if i - 1 >= 0:
        dp[i][j] += problem15(i - 1,j)
    
    if j - 1 >= 0:
        dp[i][j] += problem15(i,j - 1)
        
    return dp[i][j]

#print problem15(20,20)

def logPow(a,b):
    ret = 1
    while b > 0:
        if b & 1:
            ret *= a
        a = a*a
        b >>= 1
        
    return ret

def sumOfDigits(val):
    ret = 0
    zero = ord('0')
    s = str(val)
    for d in s:
        ret += (ord(d) - zero)
        
    return ret


def problem16():
    return sumOfDigits(logPow(2,1000))

digits = ['one','two','three','four','five','six','seven','eight','nine']
tenToTwenty = ['eleven','twelve','thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen']
timesTen = ['ten','twenty','thirty','forty','fifty','sixty','seventy','eighty','ninety']


def numberToNumberName(num):
    if num == 1000:
        return 'one thousand'
    
    if num%100 == 0:
        return digits[num/100 - 1] + ' hundred'
    
    if num > 99:
        return digits[num/100 - 1] + ' hundred and ' + numberToNumberName(num%100)
    
    if num > 9:
        if num%10 == 0:
            return timesTen[num/10 - 1]
        
        if num < 20:
            return tenToTwenty[num - 11]
        
        return timesTen[num/10 - 1] + '-' + digits[num%10 - 1]
        

    return digits[num - 1]
        

def countLetters(s):
    ret = 0
    for c in s:
        if ('a' <= c and c <= 'z'):
            ret += 1
    
    return ret   

def problem17():
    ret = 0 
    for i in range(1,1001):
        #print '{0} : {1}'.format(i, numberToNumberName(i) )
        ret += countLetters( numberToNumberName(i) )    
    
    return ret

#print problem17()

#print problem14()

#also solved problem 67
def problem18():
    fname = "input.txt"
    with open(fname) as f:
        lines = f.readlines()
        
    ret = 0
    a = []
    for line in lines:
        line = line.replace("\n","")
        a.append([int(x) for x in line.split(' ') ])

    n = len(a)
    bst = []
    bst.append([a[0][0] ])
    
    for i in range(1,n):
        bst.append([0  for i in range(0,len(a[i]))])
        m = len(a[i])
        for j in range(0,m):
            if j - 1 >= 0:
                bst[i][j] = a[i][j] + bst[i - 1][j - 1]
            if j != m - 1:
                bst[i][j] = max(bst[i][j],bst[i - 1][j] + a[i][j])
    
            
            ret = max(ret,bst[i][j])
    
    return ret
   
    
#print problem18()

def isLeapYear(year):
    if year%400 == 0:
        return True
    if year%100 == 0:
        return False
    if year%4 == 0:
        return True
    
    return False

#1 Jan 1900 was a Monday.
#How many Sundays fell on the first of the month during the twentieth century (1 Jan 1901 to 31 Dec 2000)?
def problem19():
    ret = 0
    days = ['Monday','Tuesday','Wednesday','Thursday','Friday','Saturday','Sunday']
    months = ['January','February','March','April','May','June','July','August','September','October','November','December']
    regularYearDays = [31,28,31,30,31,30,31,31,30,31,30,31]
    currentDay = 0
    for year in range(1900,2001):
        yearDays = regularYearDays[:]
        if isLeapYear(year) == True:
            yearDays[1] += 1
        
        #print 'Jan - {0} - {1}'.format(days[currentDay] , year)

        for month in range(0,12):
            #currentDay = (currentDay + yearDays[month])%7
            if days[currentDay] == 'Sunday':
                ret += 1
                print '{0}-{1}-{2}-{3}'.format(1,days[currentDay ],months[month],year)
            #for day in range(0,yearDays[month]):
                #print '{0} - {1} - {2} - {3}'.format(day + 1,days[(currentDay + day)%7],months[month],year)
            currentDay = (currentDay + yearDays[month])%7
            
    return ret

#print problem19()

def problem31():
    ways = [0 for i in range(0,201)]
    coinValues = [1,2,5,10,20,50,100,200]
    target = 200
    ways[0] = 1
    for i in range(0,8):
        j = coinValues[i]
        while j <= target:
            ways[j] += ways[j - coinValues[i]]
            j += 1
            
    return ways[target]

#print problem31()

def generateBrute(n,p,k):
    if n == 10:
        if k == 10:
            print p
        
        return
    if k > 10:
        return
    
    for i in range(0,10):
        generateBrute(n + 1,p + str(i),k)
        k += 1
    
 
def problem24():
    generateBrute(0,"",0)

def nameScore(k,name):
    ret = 0
    AIndex = ord('A')
    for c in name:
        ret += (ord(c) - AIndex + 1)
      
    ret *= k
    
    return ret

def problem22():
    fname = "input.txt"
    ret = 0
    with open(fname) as f:
        lines = f.readlines()
    
    a = []
    for line in lines:
        line = line.replace("\n","")
        line = line.replace('"',"")

        a.append([str(x) for x in line.split(',') ])    
  
    a[0].sort()
    k = 0
    for name in a[0]:
        k += 1
        ret += nameScore(k,name)
        
    return ret
    
#print problem22()

def factorial(n):
    ret = 1
    for i in range(2,n + 1):
        ret *= i
    return ret

def problem20():
    return sumOfDigits(factorial(100))



#print problem20()

D = [-1 for i in range(0,100002)]

def sumOfProperDivisors(n):
    if D[n] != -1:
        return D[n]
    
    D[n] = 0;
    
    for i in range(1,n/2 + 1):
        if n%i == 0:
            D[n] += i
        
    return D[n]

def problem21():
    ret = 0
    for i in range(1,10000):
            b = sumOfProperDivisors(i)
            if i != b and sumOfProperDivisors(b) == i:
                ret += i
    
    return ret
    
#print problem21()

def distinctPrimeFactors(n):
    ret = 0
    d = 2
    if n%d == 0:
        ret += 1
        while n%d == 0:
            n /= d
    d = 3
    
    while d*d <= n:
        if n%d == 0:
            ret += 1
            while n%d == 0:
                n /= d
        d += 2
    
    if n > 1:
        ret += 1
    
    return ret
   
def problem47():
    n = 2
    while True:
        if ( (distinctPrimeFactors(n) ^ 4) | (distinctPrimeFactors(n + 1) ^ 4) | (distinctPrimeFactors(n + 2) ^ 4) | (distinctPrimeFactors(n + 3) ^ 4)) == 0:
            return n
        n += 1
        
    return -1

#print problem47()

def isPermutationOf(a,b):
    s1 = str(a)
    s2 = str(b)
    sorted(s1, lambda x,y: cmp(x.lower(), y.lower()) or cmp(x,y))
    sorted(s2, lambda x,y: cmp(x.lower(), y.lower()) or cmp(x,y))
    return 1 if s1 == s2 else 0
    

def problem49():
    for i in range(1000,10000):
        if prime(i) == 1:
            j = 1
            while i + 2*j < 10000:
                b = i + j
                c = i + j + j
                if (isPermutationOf(i,b) and isPermutationOf(i,c)) == 1:
                    print str(i) + str(b) + str(c)
                j += 1
    

#problem49()

def problem35():
    ret = 0
    p = [0 for i in range(0,1000000)]
    for i in range(2,1000000):
        if p[i] == 0:
            j = i + i
            while j < 1000000:
                p[j] = 1
                j += i
    
    for i in range(2,100):
        n = str(i)
        l = len(n)
        ok = True
        for j in range(0,l):
            val = 0
            for k in range(0,l):
                val = val*10 + (ord(n[(j + k)%l]) - ord('0'))
            
            if p[val] == 1:
                ok = False
                break
        
        ret += ok
    return ret    
    
#print problem35()

#to do
def problem34():
    ret = 0
    fact = [0 for i in range(0,10)]
    fact[0] = 1
    for i in range(1,10):
        fact[i] = fact[i - 1]*i
        
    for i in range(10,100000000):
        factorialSum = 0
        val = i
        while val > 0:
            factorialSum += fact[val%10]
            val /= 10
            
        if i == factorialSum:
            ret += 1
    return ret

#print problem34()

def problem52():
    x = 1
    while True:
        aux = x
        digitCount = [0 for i in range(0,10)]
        while aux > 0:
            digitCount[aux % 10] += 1
            aux /= 10
        
        for i in range(2,7):
            aux = x*i
            auxCount = [0 for j in range(0,10)]
            while aux > 0:
                auxCount[aux % 10] += 1
                aux /= 10
            if auxCount != digitCount:
                break
            
        if auxCount == digitCount:
            return x
            
        x += 1
    return -1

#print problem52()

def isAbundant(x):
    s = 0
    for i in range(1,x/2 + 1):
        if x%i == 0:
            s += i
            
    return x < s

def problem23():
    a = []
    sums = [0 for i in range(0,28124)]
    for i in range(1,28124):
        if isAbundant(i):
            a.append(i)
           
    for i in range(0,len(a)):
        for j in range(i,len(a)):
            if a[i] + a[j] < 28124:
                sums[a[i] + a[j]] = 1
                
    ret = 0
    for i in range(1,28124):
        if sums[i] == 0:
            ret += i

    return ret

#print problem23()

def matrixMultiply(a,b):
    c = [[0,0],[0,0]]

    for i in range(0,2):
        for j in range(0,2):
            for k in range(0,2):
                c[i][j] += a[i][k]*b[k][j]
           
    return c

#Time complexity: O(log(n))
#fib[1] = fib[2] = 1
def fibo(n):
    n -= 2
    m = [[0,1],[1,1]]
    ret = m[:]
    while n > 0:
        if n & 1:
            ret = matrixMultiply(ret,m)
         
        m = matrixMultiply(m,m)
        n >>= 1
    
    return ret[1][1]

#Time complexity: O(log(n)*Log(n))
def problem25():
    low = 0
    high = 10000
    while low <= high:
        mid = low + (high - low)/2
        digitCount = len(str(fibo(mid)))
        if digitCount < 1000:
            low = mid + 1
        else:
            high = mid - 1
            
    return low


#print problem25()
#primitive root example --->
#The number 3 is a primitive root modulo 7 because
#we see that the period of 3^k (for k(1,6)) modulo 7 is 6. The remainders in the period,
# which are 3, 2, 6, 4, 5, 1, 
#form a rearrangement of all nonzero remainders modulo 7, implying that 3 is indeed a primitive root modulo 7

#If 10 is a primitive root modulo p, the period is equal to p - 1
#; if not, the period is a factor of p - 1
#This result can be deduced from Fermat's little theorem, which states that 10^(p - 1) = 1 (mod p)

def problem26():
    ret = 0
    for i in range(1,1000):
        if prime(i):
            use = [0 for j in range(0,i)]
            power10 = 1
            use[0] = 1
            for j in range(1,i):
                power10 = power10*10%i
                use[power10] = 1
            if use.count(1) == i:
                ret = max(ret,i - 1)
                #print i
    return ret

#print problem26()

def problem53():
    ret = 0
    C = [[1,0]]
    for i in range(1,101):
        C.append([0 for k in range(0,i + 2)])
        C[i][0] = 1
        for j in range(1,i + 1):
            C[i][j] = C[i - 1][j] + C[i - 1][j - 1]
            ret += (C[i][j] > 1000000)
     
    #print C       
    return ret      

#print problem53()

def problem81():
    n = m = 80
    fname = "input.txt"
    with open(fname) as f:
        lines = f.readlines()
        
    a = []
    for line in lines:
        line = line.replace('\n','')
     
        a.append([int(x) for x in line.split(',') ])
        
    bst = [[10000000000 for j in range(0,m)] for i in range(0,n)]
   
    for i in range(0,n):
        for j in range(0,m):
            if i == 0 and j == 0:
                bst[i][j] = a[i][j]
                continue
            if i == 0:
                bst[i][j] = bst[i][j - 1] + a[i][j]
            else: 
                if j == 0:
                    bst[i][j] = bst[i - 1][j] + a[i][j]
                else:
                    bst[i][j] = min(bst[i - 1][j],bst[i][j - 1]) + a[i][j]
         
    
    return bst[n - 1][m - 1]           
    
#print problem81()
 
def problem82():
    n = m = 80
    fname = "input.txt"
    with open(fname) as f:
        lines = f.readlines()
        
    a = []
    for line in lines:
        line = line.replace('\n','')
     
        a.append([int(x) for x in line.split(',') ])
        
    bst = [[10000000000 for j in range(0,m)] for i in range(0,n)]
    for j in range(0,m):
        if j == 0:
            for i in range(0,n):
                bst[i][j] = a[i][j]
        else:
            bst[0][j] = bst[0][j - 1] + a[0][j]
            for i in range(1,n):
                bst[i][j] = min(bst[i - 1][j],bst[i][j - 1]) + a[i][j]
            
            for i in range(1,n):
                ri = n - i - 1
                bst[ri][j] = min(bst[ri][j],bst[ri + 1][j] + a[ri][j])
      
        
    ret = 100000000000
    for i in range(0,n):
        ret = min(ret,bst[i][m - 1])
    
    return ret

#print problem82()

def inRange(w):
    return ( 0 <= w[0] and w[0] < 80 and 0 <= w[1] and w[1] < 80)
 
def problem83():               
    n = m = 80
    dx = [-1,1,0,0]
    dy = [0,0,-1,1]
    fname = "input.txt"
    with open(fname) as f:
        lines = f.readlines()
        
    a = []
    for line in lines:
        line = line.replace('\n','')
     
        a.append([int(x) for x in line.split(',') ])
        
    bst = [[10000000000 for j in range(0,m)] for i in range(0,n)]
    
    q = deque()
    q.append([0,0])
    bst[0][0] = a[0][0]
    while q:
        v = q.popleft()
        for k in range(0,4):
            w = [v[0] + dx[k],v[1] + dy[k]]
            if inRange(w) and bst[w[0]][w[1]] > bst[v[0]][v[1]] + a[w[0]][w[1]]:
                bst[w[0]][w[1]] = bst[v[0]][v[1]] + a[w[0]][w[1]]
                q.append(w)
        
   
    return bst[n - 1][m - 1] 
#print problem83()

def problem29():
    s = set()
    for a in range(2,101):
        for b in range(2,101):
            s.add(logPow(a,b))
    
    return len(s)
    
#print problem29()

def problem48():
    sum = 0
    mod = 10000000000
    for i in range(1,1001):
        sum += logPow(i,i)
        sum %= mod
    
    return sum

#print problem48()


def problem49_():
    for i in range(1000,10000):
        if i == 1487:
            continue
        j = 1;
        while  i + (j<<1) < 10000:
            if isPermutationOf(i,i + (j<<1)) == 1:
                return i
            j += 1
    return -1;    

def problem51():
    n = 1000000
    p = [0 for i in range(n)]
    i = 3
    while i*i < n:
        if p[i] == 0: 
            j = i + i + i
            while j < n:
                p[j] = 1
                j += i<<1
                
        i += 2         
        
    i = 3
    primes = set()
    sums = []
    sums.append(0)
    sums.append(2)
    primes.add(2)
    k = 1
    while i < n:
        if p[i] == 0:
            primes.add(i)
            k += 1
            sums.append(i)
            sums[k] += sums[k - 1]
        i += 2   
        
    dim = len(sums)
    k = 1
    numPrimes = 0
    primeNum = -1
    while k < dim:
        i = k
        while i > 0 and sums[k] - sums[i - 1] < n:
            if sums[k] - sums[i - 1] in primes and k - i + 1 > numPrimes:
                numPrimes = k - i + 1
                primeNum = sums[k] - sums[i - 1]
            i -= 1
        i = k
        while i < dim and sums[i] - sums[k - 1] < n:
            if sums[k] - sums[i - 1] in primes and i - k + 1 > numPrimes:
                numPrimes = i - k + 1
                primeNum = sums[i] - sums[k - 1]
            i += 1
        k += 1     
    print '{0} {1}'.format(primeNum,numPrimes)
    
    
#problem51()

def superx():
    
    for i in range(100):
        print (i)
                



class Shape:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    description = "This shape has not been described yet"
    author = "Nobody has claimed to make this shape yet"
    def area(self):
        return self.x * self.y
    def perimeter(self):
        return 2 * self.x + 2 * self.y
    def describe(self,text):
        self.description = text
    def authorName(self,text):
        self.author = text
    def scaleSize(self,scale):
        self.x = self.x * scale
        self.y = self.y * scale
    

rectangle = Shape(100, 45)

#print rectangle.area()


    