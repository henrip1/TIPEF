from scipy.optimize import minimize
import l_json

def calculate_error(params, S0, E0, I0, R10, R20, donnes):
    a, b, c, f = params
    _, _, I, _, _ = equations_SEIR(y, days, a, b, c, f)
<<<<<<< HEAD
    IR,_,_ = donnes
    error = np.mean((I-IR)**2)
=======
    error = np.mean((S-external_date)**2) + np.mean((I-external_date)**2)
>>>>>>> d8fd8a4a48567ad0871a6176c963c2e987904723
    return error

S0 = 999999
E0 = 1
I0 = 0
R10 = 0
R20 = 0
y0 = S0, E0, I0, R10, R20
a_v = []
b_v = []
c_v = []
f_v = []
day = 200
ini_params = [0.3,0.1,0.2,0.8]

for day in range(1, day + 1):

    result = minimize(calculate_error, ini_params, args=( S0, E0, I0, R10, R20, day, donnes_nuplets[day+1] ))#jour un par

    a_opt, b_opt, c_opt, f_opt = result.x

    a_v.append(a_opt)
    b_v.append(b_opt)
    c_v.append(c_opt)
    f_v.append(f_opt)



    S, R, I, R1, R2 = SEIR_model(y, 2, a_opt, b_opt, c_opt, f_opt)
    
        
    ini_params = [a_opt, b_opt, c_opt, f_opt]
    
    plt.figure(figsize=(12,8))
    
    plt.subplot(2,2,1)
    plt.plot(S_prediction, label='S (modèle SEIR)')
    plt.plot(external_date, label='Données réelles (S)', linestyle='dashed')
    plt.xlabel('Jours')
    plt.ylabel('Nombre de personnes')
    plt.title('Sains (S)')
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(S_prediction, label='E (modèle SEIR)')
    plt.plot(external_date, label='Données réelles (E)', linestyle='dashed')
    plt.xlabel('Jours')
    plt.ylabel('Nombre de personnes')
    plt.title('Exposés (E)')
    plt.legend()

    plt.tight_layout()
    plt.show()

