import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def système(t, Y):

    T = Y[0]
    I = Y[1]

    dTdt = ((RL + Rtes(T, I))*I**2 - Pbain(T) + Pray(t))/C
    dIdt = (V - I*(RL + Rtes(T, I)))/L
    return [dTdt, dIdt]


def get_Rtes(alpha,beta,I0,T0,R0):
  def Rtes(T,I):
    return R0 + alpha*(R0/T0)*(T-T0) + beta*(R0/I0)*(I-I0)
  return Rtes

def get_Pbain(G, Tbain):
  def Pbain(T):
    return G*(T-Tbain)
  return Pbain

def bruit():
  while True:
    sample = np.random.randn()  
    noise = sample / 2 
    if noise >= -0.4 and noise <= 0.4:  
        return 1 + noise

def escalier(t):
  if t <= 2:
    Pray = 1
  elif 2 < t < 2.2:
    Pray = 1.25
  elif 2.6 <= t < 2.9:
    Pray = 0.75
  elif 2.9 <= t < 3.2:
    Pray = 0.9
  elif 5 <= t < 5.6:
    Pray = 1.5
  elif 5.6 <= t < 6.5:
    Pray = 1.3 
  elif 7 <= t < 7.5:
    Pray = 0.5
  else :
    Pray = 1
  return Pray

def pt_isolés(t):
  if np.isclose(t, 1, atol=0.001) :
    Pray = 2
  elif np.isclose(t, 3, atol=0.001):
    Pray = 1
  elif np.isclose(t, 6, atol=0.001):
    Pray = 1.5
  else :
    Pray = 0
  return Pray

def bruit_poisson(x):
  attenuation = 0.2
  noise = (np.random.poisson(lam = x) - x)* attenuation
  return noise


def Pray(t):
    #Pray = 0 if t < 5 else 2
    #Pray = 1
    #Pray = 2 if t < 5 else 0
    #Pray = np.exp(-1*t)
    #Pray =np.sin(t) + bruit()
    Pray = 2 if 4 < t < 6 else 0
    #Pray = 2 if np.isclose(t, 1, atol=0.001) else 0

    #Pray = pt_isolés(t)
  
    #Pray = bruit()

    #Pray = 2 + bruit() if 4 < t < 6 else 0 + bruit()
    #Pray = escalier(t)
    #Pray1 = escalier(t)
    #Pray = Pray1 + bruit_poisson(Pray1)

    #Pray1 = 2 if 4 < t < 6 else 0
    #Pray = Pray1 + bruit_poisson(Pray1)

    #Pray = 0 if t < 5 else 1
    #Pray = Pray1 + bruit_poisson(Pray1)
  
    return Pray

t_span = (0, 10)
T0 = 2
I0 = 1
R0 = 0.1
alpha = 1
beta = 1
Tbain = 1
G = 1
C = 1
V = 1
RL = 0.9
L = 1

T_0 = 2
I_0 = 1



Rtes = get_Rtes(alpha, beta, I0, T0, R0)
Pbain = get_Pbain(G, Tbain)

sol = solve_ivp(
  système, t_span, 
  [T_0, I_0],  
  dense_output=True,
  max_step=0.01,  
  first_step=0.01 
)

t = sol.t
T = sol.y[0]
I = sol.y[1]

intervalle_I = np.array([])
intervalle_P = np.array([])
index_t_5 = np.where(np.isclose(t, 5, atol=0.001))[0]

deltaI = abs(I[-1] - I[index_t_5])

Rtes_values = get_Rtes(alpha, beta, I0, T0, R0)(T, I)
Pray_values = [Pray(ti) for ti in t]

max_Rtes = max(Rtes_values)
max_Rtes_index = np.argmax(Rtes_values)
t_max_Rtes = t[max_Rtes_index]
t_max_Rtes_arrondi = np.ceil(t_max_Rtes * 100) / 100

min_I = min(I)
min_I_index = np.argmin(I)
t_min_I = t[min_I_index]
t_min_I_arrondi = np.ceil(t_min_I * 100) / 100

max_T = max(T)
max_T_index = np.argmax(T)
t_max_T = t[max_T_index]
t_max_T_arrondi = np.ceil(t_max_T * 100) / 100

delta_I = I - I0
delta_R = Rtes_values - R0

print('Rtes final',Rtes_values[-1])
print('I final',I[-1])
print('T final',T[-1])
print('La valeur max de Rtes est',max_Rtes, 'pour t =', t_max_Rtes_arrondi)
print('La valeur min de I est',min_I, 'pour t =',t_min_I_arrondi)
print('La valeur max de T est',max_T, 'pour t =', t_max_T_arrondi)


plt.subplot(3, 1, 1)
plt.plot(t, T, label="T en K", color = 'blue')
plt.plot(t, I, label="I en A", color='orange')
plt.plot(t, Rtes_values, label="$R_{TES}$ en $\Omega$", color ='red')
plt.ylabel("T(t), I(t) et $R(t)$")
plt.ylabel('I(t) en A')
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(t, Pray_values, label="$P_{ray}$ en W", color="red")
plt.ylabel("Pray(t) en W")
plt.xlabel('t en s')
plt.legend()

plt.subplot(3, 1, 3)
plt.plot(t, delta_I, label='$\delta I$ en A', color ='orange')
plt.plot(t, delta_R, label='$\delta R$ en $\Omega$', color='red')
plt.axhline(y=0, color='gray', linestyle='--')
plt.ylabel('$\delta I$ et $\delta R$')
plt.xlabel('t en s')
plt.legend()
plt.savefig("equadiff_projet_normalisé.png")
plt.close()


