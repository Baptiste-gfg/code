import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import ScalarFormatter
from scipy.interpolate import interp1d

def système(t, Y):
    T = Y[0]
    I = Y[1]
    dTdt = ((RL + Rtes(T, I))*I**2 - Pbain(T) + Pray(t))/C
    dIdt = (V - I*(RL + Rtes(T, I)))/L
    return [dTdt, dIdt]

def get_Rtes(alpha, beta, I0, T0, R0):
    def Rtes(T, I):
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

def bruit_poisson(x):
    attenuation = 0.001
    noise = (np.random.poisson(lam = x) - x)* attenuation
    return noise

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
intervalle_I = np.array([])
intervalle_P = np.array([])
t = np.array([])

start = 1  
stop = 100  
num_points = 200  
deltaP_values = np.logspace(-6,-2, num=num_points)

intervalle_I = np.array([])
intervalle_P = np.array([])

for Delta_P in deltaP_values:
    def Pray(t):
        Pray1 = 1 if t < 5 else 1 + Delta_P
        Pray = Pray1 + bruit_poisson(Pray1)
        return Pray

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
    index_t_5 = np.where(np.isclose(t, 5, atol=0.001))[0]
    index_t_car = np.where(np.isclose(t, 7.8, atol=0.001))[0]
    
    Delta_I = abs(I[index_t_car] - I[index_t_5])

    intervalle_I = np.append(intervalle_I, Delta_I)
    intervalle_P = np.append(intervalle_P, Delta_P)

interpolateur = interp1d(intervalle_P, intervalle_I, kind='linear', fill_value='extrapolate')
delta_P_target = 1e-4
delta_I_estime = interpolateur(delta_P_target)


first_value = intervalle_I[0]
ecart = abs((delta_I_estime - first_value)/ first_value)
print(f"Pour ΔP = {delta_P_target:.1e}, ΔI = {delta_I_estime:.5e}")

def format_power(x, pos):
    return f'{x:.2e}'

print('l\'écart relatif vaut',ecart)
plt.loglog(intervalle_P, intervalle_I, label='Données')
plt.axhline(y=first_value, color='red', linestyle='--', label='Asymptote')
plt.xlabel(r'$\Delta P$ en W')
plt.ylabel(r'$\Delta I$ en A')
plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
plt.tight_layout()
plt.grid(True)
plt.legend()
plt.savefig("equadiff_projet_normalisé.png")
plt.close()