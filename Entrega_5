import scipy as sp

def ynormal(B,ss,yini,tol,C0,S0,n):
    yn = yini
    q=30 #m**3/s
    # print "Parto en: {}".format(yn)

    while tol>(10**(-6)):
        #Calculamos Variables
        A = B*yn + ss*(yn**2)
        P = B + 2*yn*sp.sqrt(ss**2+1)
        dAdyn = B + 2*ss*yn
        dPdyn = 2*sp.sqrt(ss**2+1)
        Q = (C0*(A**(5./3.)*(P**(-2./3.))*sp.sqrt(S0))/n)

        #Preparara Newton-Raphson
        fyn = (A**(5./3.))*(P**(-2./3.))-(n*q)/(C0*sp.sqrt(S0))
        dfyn = (5./3.)*(A**(2./3.))*(P**(-2./3.))*dAdyn - (2./3.)*(A**(5./3.))*(P**(-5./3.))*dPdyn
        yn2 = yn - (fyn/dfyn)


        #Calculo tolerancia
        tol = abs(yn2-yn)

        #Se actualiza yn
        yn = yn2
        # print "voy en: {}".format(yn)
    return yn

def NRalfa(Q,B,ss,n,Co,S):
	#escrbir esta funcion segun el video


N =1000       #numero de nodos
L = 150000    #Longitud del rio 
B = 100       #Ancho base canal
S = 0.001     #pendiente del canal
n = 0.045     # coef manning
ss = 0	      # Side slope
C0 = 1.485    # unidades british
NC = 1.	      # numero de curant
Tfin = 600*60 # tiempo total de simulacion en segundos
g = 32.2      # gravedad

dx = L/(N-1)
x = sp.zeros(N)

for i in range(len(x)):
	x[i]=(i-1)*dx 


#Condiciones iniciales k=1 t=0
Q = sp.zeros((N,N))
Q[1,:]=250.

y = sp.zeros((N,N))
y[1,:] = ynormal(1,B,Q[1,1],ss,n,S,C0)

A = sp.zeros((N,N))
A[1,:] = B*y[1,1]

V = sp.zeros((N,N))
V[1,:] = Q[1,1]/A[1,1]

#simulacion

t = 0
k = 0
pi = sp.pi 

while t<Tfin:
	
	dt=NC*dx/(mean(V[k,:]));
	k+=1
	t+=dt
	tshow = t/60

	if t<=150*60:
		Q[k,1] = 250+750/pi*(1-cos(pi*t/(60*75)))
	else:
		Q[k,1]=250

	ay = NRalfa(Q[k,1],B,ss,n,Co,S)  #NR  para  encontar  valor  de  A  e  Y  iterando
	A[k,1] = ay[0]
	y[k,1] = ay[1]
	V[k,1]  =  Q[k,1]/A[k,1]

	# Moviendo  estencil  en  el  tiempo  k,  desde  nodo  1  al  nodo  N
	for i in range(len(x)-1):

		Q[k,i+1]=Q[k,i]-dt/dx*(A([k,i]-A[k-1,i])    #ecuacion  de  continuidad
		ay = NRalfa(Q[k,i+1],B,ss,n,Co,S)
		A[k,i+1] = ay[0]
		y[k,i+1] = ay[1]
		V[k,i+1] = Q[k,i+1]/A[k,i+1]















