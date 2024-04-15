import numpy as np

def u(x, y):
    return x*(x-1)*y*(y-1)
    #return x*(x-2)*y*(y-2)

def f(vn,i,j):
    return A(vn)[i][j] - b(i*h1,j*h2)
    #чтобы проверить ур пуассона в прямоуг - надо возвращать в f то что снизу(f= -лапласU)
    #return -2*(x**2-x)-2*(y**2-y)

    #return -2*(x**2-2*x)-2*(y**2-2*y)

def mu_k2(vn, i, k2):
    global h1, h2, l1, l2, N1, N2
    sum=0
    for j in range(1,N2):
        sum+=f(vn, i, j)*np.sin(k2*np.pi*j/N2)
        
    return sum
def mu_k1k2(vn, k1, k2):
    global h1, h2, l1, l2, N1, N2
    sum=0
    for i in range(1,N1):
        sum+=mu_k2(vn, i, k2)*np.sin(k1*np.pi*i/N1)
    return sum
def v_k2(vn, i, k2):
    global h1, h2, l1, l2, N1, N2
    sum=0
    for k1 in range(1,N1):
        sum+= mu_k1k2(vn, k1, k2) / sobstv_ch(k1, k2)*np.sin(k1*np.pi*i/N1)
    return sum
def sobstv_ch(k1, k2):
    global h1, h2, l1, l2, N1, N2
    return (4./h1**2)*np.sin(k1*np.pi*h1/(2.*l1))**2 + (4./h2**2)*np.sin(k2*np.pi*h2/(2*l2))**2 
def v(vn, i,j):
    global h1, h2, l1, l2, N1, N2
    sum=0

    for k2 in range(1,N2):
        sum+=v_k2(vn, i, k2)*np.sin(k2*np.pi*j/N2)
    return (4/(N1*N2))*sum

def k11(x,y):
    return 1+x**2+y**2
def k22():
    return 1

#lambda_xy=lambda_yx=0
#k11(x,y) k22=1 k12=k21=0

def lambda_xx(vn, i, j):
    #if j-1 == 0:
    #    vn[i+1][j-1] = 0
    #    vn[i][j-1] = 0
    #elif j + 1 == N1:
    #    vn[i][j+1] = 0
    #    vn[i-1][j+1] = 0
    #elif i - 1 == 0:
    #    vn[i-1][j+1] = 0
    #    vn[i-1][j] = 0
    #elif i + 1 == N1:
    #    vn[i+1][j] = 0
    #    vn[i+1][j-1] = 0
    return (1/h1**2)*(vn[i+1][j]*0.5*(k11(i*h1, j*h2) + k11((i+1)*h1, j*h2)) - vn[i][j]*0.5*( 2*k11(i*h1, j*h2) + k11((i+1)*h1, j*h2)+ k11((i-1)*h1, j*h2) ) + vn[i-1][j]*0.5*(k11(i*h1, j*h2) + k11((i-1)*h1, j*h2)))

def lambda_yy(vn, i, j):
    #if j-1 == 0:
    #    vn[i+1][j-1] = 0
    #    vn[i][j-1] = 0
    #elif j + 1 == N1:
    #    vn[i][j+1] = 0
    #    vn[i-1][j+1] = 0
    #elif i - 1 == 0:
    #    vn[i-1][j+1] = 0
    #    vn[i-1][j] = 0
    #elif i + 1 == N1:
    #    vn[i+1][j] = 0
    #    vn[i+1][j-1] = 0
    return (1/h2**2)*(vn[i][j+1] - vn[i][j]*2. + vn[i][j-1])

#правая часть урния Au=f с K11=1+x^2+y^2 K22=1 u=x(x-1)y(y-1)
def b(x,y):
    return -2*((y**2-y)*(3*x**2+y**2-x+1)+ x**2-x)

def A(vn):
	A_vn = np.zeros((N1+1,N2+1))
	for i in range(1, N1):
		for j in range(1, N2):
			A_vn[i][j] = -(lambda_xx(vn, i, j) + lambda_yy(vn, i, j) )
	return A_vn

def alpha():
    #ищем через M m - ограничения для K11
    #для k11=1+x**2+y**2 на ед квадрате:
    c1=1.
    c2=3.
    return (c2-c1)/(c2+c1)


#N колво h внутри отрезка(потому делим l рочно на N), те N+1 точка будет 
N1=5
N2=5
#квадрат [0 l1] [0 l2]
l1=1. #2
l2=1. #2
h1=l1/(N1)
h2=l2/(N2)

eps=10**(-3)

v0 = np.zeros((N1+1,N2+1))
vn = v0

x=np.linspace(0,l1,N1+1)
y=np.linspace(0,l2,N2+1)
X,Y=np.meshgrid(x,y)
right_p = b(X, Y)

#решение с кт будем сравнивать получившуюся vn
etalon = u(X, Y)

while(np.linalg.norm(A(vn)-right_p, ord=np.inf ) > eps ):
    # решаем By=A(vn)-b 
    y = np.zeros((N1+1,N2+1)) #от 0 до N1 индексы
    for i in range(1, N1):
        for j in range(1, N2):
            y[i][j]=v(vn, i, j)

    print(f'norma rn:{np.linalg.norm(A(vn)-right_p, ord=np.inf )} ')
    #print(f'||vn-u||={np.linalg.norm(vn-etalon, ord=np.inf )}')

    #попытка написать вчера вечером поиск альфа третьим методом, как у полины - мб есть ошибки но не должно быть:
    #print(f'vn до:{vn}')
    ##
    #rnn=A(vn)-right_p 
    #sum1=0.
    #sum2=0.
    #for i in range(1,N1):
    #    for j in range(1,N2):
    #        sum1+=rnn[i][j]*y[i][j]*h1**2

    #for i in range(1,N1):
    #    for j in range(1,N2):
    #        sum2+=A(y)[i][j]*y[i][j]*h1**2

    #alpha=sum1/sum2*10**(2)
    ##
    vn -= alpha()*y
    #vn -= alpha*y


print(f'u на сетке:\n {etalon}')
print(f'vn :\n {vn}')

print(f'||vn-u||={np.linalg.norm(vn-etalon, ord=np.inf )}')



