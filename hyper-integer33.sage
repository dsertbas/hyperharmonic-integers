# The following SageMath code is designed to solve this specific question: "Are there any hyperharmonic integers?" 
# Therefore, it is not optimized and may contain some flaws or bugs in it. 
# Nevertheless, it is sufficient to obtain the results that are described in D. C. Sertbas, "Hyperharmonic integers exist", submitted.

# r=64*(2^alpha-1)+32 --> Main Aim: Finding appropriate alpha values for all primes p<=33.

#  ------------------  The cases when p in {7,13,17,19,23,29,31}  ------------------ #

def pos_rps(p):		#computes possible r (mod p)'s for the primes bigger than floor(n/p), where p does not divide 33
	if floor(33/p)>=p or 33%p==0 or not is_prime(p):   
		return []
	n_p=Mod(33,p)
	R=[1-b for b in [1..n_p]]
	R.sort()
	return R

def dlog_two(n,p):		#computes the corresponding discrete logarithm in base 2. We will get a solution in mod(ord_2(p))
# p is not necessarily a prime number, as we will see in the case when p=3 (there p=3^6).
	if gcd(2,p)==2 or gcd(n,p)!=1:
		return p	
# This result is incorrect. But we want to have some numerical results for these cases (the cases which do not return any solution), in order to eliminate them easily later. (*)
# So for these cases, we assign a value which does not lie in the range of dlog_2 (discrete log in base 2)
	elif primitive_root(p, check=False)==2:
		return discrete_log(Mod(n,p),Mod(2,p))
	else:
		ord=Mod(2,p).multiplicative_order()
		powers=[Mod(2^a,p) for a in range(ord)]
		if not n in powers:		# The same reasoning in (*) occurs also here.
			return p		# So, we assigned the value p again.
		return powers.index(n)


P_first=[7,13,17,19,23,29,31]

def pos_alphaps(p):	#finding all possible alphas with respect to prime p in {7,13,17,19,23,29,31}
	if not p in P_first:   
		return False
	k=Mod(2^(-6),p)
	A=[Mod(x,p) for x in pos_rps(p)]			# the list of 1-b_p's in (mod p)
	D=[Mod(1+k*(a-32),p) for a in A]			# therefore, this part looks different from the paper, but the resulting list of alphas are the same 
	S=[dlog_two(d,p) for d in D if dlog_two(d,p)!=p]
	al=set(S)						# unify all overlapping solutions
	AL=list(al)
	AL.sort()
	return AL


first_alphas=[pos_alphaps(p) for p in P_first]
first_ords=[Mod(2,p).multiplicative_order() for p in P_first]		#the list of orders of 2 in (ZZ/pZZ)^\times for p in {7,13,17,19,23,29,31}

def common_crt(Q,W):	
# Finding all common solutions (mod lcm(W)) where each Q[i] is a list of possible values in (mod W[i]).
# SageMath does not contain such a built-in function. So, it was needed to be written.
	l=len(W)
	T=[]
	parameter=0
	for q in Q:
		if type(q) is not list:
			parameter=1
	if l==0 or len(Q)!=l or [w for w in W if w==1]!=[] or parameter==1: 
# all unwanted cases are eliminated
		return [T,1]
	elif l==1:
		return [Q,lcm(W)]
	elif l==2:
		if Q[0]==[] or Q[1]==[]:
			return [T,lcm(W)]
		else:
			d=gcd(W[0],W[1])
			for x in Q[0]:
				for y in Q[1]:
					if (x-y)%d==0:
						t=crt(x,y,W[0],W[1])
						T.append(t)
			T.sort()
			return [T,lcm(W)]
	else:
		A=[Q[l-2],Q[l-1]]
		B=[W[l-2],W[l-1]]
		C=common_crt(A,B)
		Q1=[Q[n] for n in range(l-2)]
		Q1.append(C[0])
		W1=[W[n] for n in range(l-2)]
		W1.append(C[1])
		CC=common_crt(Q1,W1)
		c=set(CC[0])
		C0=list(c)
		C0.sort()
		return ([C0,CC[1]])

CAL33_first=common_crt(first_alphas,first_ords) # CAL33_first = AL33_first

# lcm(first_ords) = 27720 ( = CAL33_first[1] )

A_0 = CAL33_first[0]

# len(A_0) = 196

AL33_first=[[79,  139,  359,  699,  839,  1059,  1259,  1339,  1399,  1539,  1619,  1674,  1759,  1899,  1954,  2239,  2459,  2514,  2599,  2659,  2799,  2879,  3019,  3159,  3219,  3354,  3779,  3919,  4139,  4194,  4419,  4474,  4479,  4619,  4759,  4839,  4979,  5034,  5319,  5539,  5679,  5739,  5879,  6019,  6099,  6239,  6299,  6379,  6999,  7279,  7499,  7554,  7559,  7639,  7699,  7839,  7919,  8059,  8399,  8539,  8619,  8759,  8819,  8959,  9099,  9179,  9234,  9319,  9459,  9514,  10079,  10159,  10359,  10579,  10719,  10779,  10914,  10919,  11139,  11479,  11619,  11699,  11754,  11839,  12039,  12179,  12319,  12399,  12539,  12594,  12679,  13239,  13434,  13439,  13579,  13659,  13799,  13859,  13939,  13999,  14219,  14274,  14559,  14699,  14919,  15119,  15199,  15259,  15399,  15479,  15619,  15759,  16099,  16319,  16459,  16519,  16659,  16739,  16794,  16879,  17019,  17074,  17079,  17639,  17779,  17999,  18279,  18339,  18474,  18479,  18619,  18699,  18839,  19179,  19399,  19539,  19594,  19599,  19739,  19879,  19959,  20099,  20154,  20159,  20239,  20859,  20994,  21139,  21359,  21419,  21499,  21559,  21699,  21779,  21834,  21919,  22259,  22399,  22479,  22619,  22674,  22679,  22819,  22959,  23039,  23179,  23319,  23514,  23939,  24019,  24219,  24439,  24579,  24634,  24639,  24779,  24999,  25339,  25479,  25559,  25699,  25899,  26034,  26039,  26179,  26259,  26399,  26539,  27099,  27154,  27299,  27439,  27519,  27659,  27714,  27719], 27720]

#  ------------------  The case when p = 11  ------------------ #

R11=[1,2,3,4,5,6,7,8,9,10,11,78,79,80,81,82,83,84,85,86,87,88]

A_00=[alpha for alpha in A_0 if Mod(64*(2^alpha-1)+32,121) in R11] 
# This is the set that is given as A_0' in D. C. Sertbas, "Hyperharmonic integers exist".

# len(A_00) = 52

# A_00 = AL33_second[0] 

AL33_second = [[1259, 1339, 2659, 2799, 2879, 3019, 4194, 4419, 5739, 5879, 6099, 7279, 7499, 7639, 8819, 8959, 9179, 10359, 10579, 10719, 12039, 13434, 13439, 13579, 13659, 13799, 15119, 15199, 16519, 16659, 16739, 16879, 18279, 19594, 19599, 19739, 19959, 21139, 21359, 21499, 22674, 22679, 22819, 23039, 24219, 24439, 24579, 25899, 27299, 27439, 27519, 27659], 27720]

#  ------------------  The case when p = 5  ------------------ #

# A_00mod25=[Mod(64*(2^alpha-1)+32,25) for alpha in A_00]

bes=[a for a in A_00 if not Mod(64*(2^alpha-1)+32,25) in [0,19]] # bes=[]

def corresponding_numerator(n): 	
# finds the corresponding numerator of h_{33}^{(r)} in terms of x, when we consider the p-adic valuation of h_{33}^{(r)}, for p=3,5 (check the explanation below for the details).
	x=PolynomialRing(ZZ,'x').gen()
	f=1
	for i in range(n):
		f=f*(x+i)
	g=derivative(f,x)
	return g

# Here, we see that there are 7 multiples of 5 in {r,r+1,...,r+32}, when r=0,19 (mod 25).
# Therefore, n needs to be taken as 7 in the case of p=5.
# Similarly, there are 33/3=11 multiples of 3 in {r,r+1,...,r+32}.
# Thus, we need to consider n=11 for the latter case.

g_5 = corresponding_numerator(7)
g_3 = corresponding_numerator(11)

# g_5(x) = 7*x^6 + 126*x^5 + 875*x^4 + 2940*x^3 + 4872*x^2 + 3528*x + 720

# g_3(x) = 11*x^10 + 550*x^9 + 11880*x^8 + 145200*x^7 + 1104411*x^6 + 5412330*x^5 + 17084650*x^4 + 33638000*x^3 + 38260728*x^2 + 21257280*x + 3628800

sifir=0			#g_5(sifir)%5=0
ondokuz=19		#g_5(ondokuz)%5=0

# If r = 0 (mod 25), then g_5(c_5)=0 (mod 5), as c_5 = ceil(r/5).
# Similarly if r = 19 (mod 25), then g_5(c_5)=0 (mod 5).

#  ------------------  The case when p = 3  ------------------ #

R=QuotientRing(ZZ,3^5*ZZ)
xx=PolynomialRing(R,'xx').gen()
ff=1
for i in range(11):
	ff=ff*(xx+i)
gg=derivative(ff,xx)		# gg(x) = g_3(x) (mod 3^5) where the coefficients of gg(x) in [0..(3^5-1)]
pos_roots3=[xx for xx in range(3^5) if gg(xx)==0]

# len(pos_roots3) = 18

pos_r3s=[]
for i in range(3):
	A=[3*s-i for s in pos_roots3]
	pos_r3s.extend(A)
pos_r3s.sort()

k_3=Mod(2^(-6),3^6)
S_3=[Mod(x,3^6) for x in pos_r3s]
D_3=[Mod(1+k_3*(s_3-32),3^6) for s_3 in S_3]
A_3=[dlog_two(d,3^6) for d in D_3 if dlog_two(d,3^6)!=3^6]
al_3=set(A_3)	# unify all overlapping solutions
AL_3=list(al_3)
AL_3.sort()

# len(AL_3) = 36

al_3mod = Mod(2,3^6).multiplicative_order()

# al_3mod = 3^6-3^5 = 486, as 2 is a primitive root modulo 3^6

#  ------------------  Combining all cases  ------------------ #

ALPHAS = common_crt([AL_3,A_00],[al_3mod,AL33_second[1]]) # alphas33 = ALPHAS[0]

# Recall that A_00 = AL33_second[0] = A_0'

alphas33=[2659, 23039, 28979, 30599, 30739, 36539, 36679, 38299, 44239, 64619, 70559, 72179,  72319,  78114,  78119,  78259,  79879,  85819,  106199,  112139,  113759,  113899,  119699,  119839,  121459,  127399,  147779,  153719,  155339,  155479,  161274,  161279,  161419,  163039,  168979,  189359,  195299,  196919,  197059,  202859,  202999,  204619,  210559,  230939,  236879,  238499,  238639,  244434,  244439,  244579,  246199,  252139,  272519,  278459,  280079,  280219,  286019,  286159,  287779,  293719,  314099,  320039,  321659,  321799,  327594,  327599,  327739,  329359,  335299,  355679,  361619,  363239,  363379,  369179,  369319,  370939,  376879,  397259,  403199,  404819,  404959,  410754,  410759,  410899,  412519,  418459,  438839,  444779,  446399,  446539,  452339,  452479,  454099,  460039,  480419,  486359,  487979,  488119,  493914,  493919,  494059,  495679,  501619,  521999,  527939,  529559,  529699,  535499,  535639,  537259,  543199,  563579,  569519,  571139,  571279,  577074,  577079,  577219,  578839,  584779,  605159,  611099,  612719,  612859,  618659,  618799,  620419,  626359,  646739,  652679,  654299,  654439,  660234,  660239,  660379,  661999,  667939,  688319,  694259,  695879,  696019,  701819,  701959,  703579,  709519,  729899,  735839,  737459,  737599,  743394,  743399,  743539,  745159]

alphas33_mod = 748440 # alphas33_mod = ALPHAS[1]

# len(alphas33) = 153

#  ------------------  TESTING  ------------------ #

def hyperharmonic(n,r):			#computing h_n^{(r)}
	A=binomial(n+r-1,n)
	a=0
	for i in range(r,n+r):
		a=a+(1/i)
	return A*a 

# [alpha for alpha in alphas33 if not hyperharmonic(33,64*(2^alpha-1)+32) in ZZ] = []
