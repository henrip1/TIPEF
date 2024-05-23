
S0, E0, I0, R10, R20 = y

def calculate_error(params, S0, E0, I0, R10, R20, external_date):
    a, b, c, f = params
    _, _, I, _, _ = equations_SEIR(y, days, a, b, c, f)
    error = np.mean((S-external_date)**2) + np.mean((I-external_date)**2)
    return error

a_v = []
b_v = []
c_v = []
f_v = []


for day in range(1, days + 1):

    result = minimize(calculate_error, ini_params, args=( S0, E0, I0, R10, R20, day, external_data[]))#jour un par

    a_opt, b_opt, c_opt, f_opt = result.x

    a_v.append(a_opt)
    b_v.append(b_opt)
    c_v.append(c_opt)
    f_v.append(f_opt)



    S, R, I, R1, R2 = SEIR_model(y, 2, a_opt, b_opt, c_opt, f_opt)  
        
    ini_parems = [a_opt, b_opt, c_opt, f_opt]

    plt.figure(figsize=(12,8))
    
    plt.subplot(2,2,1)
    plt.plot(S_prediction, label='S (modèle SEIR)')
    plt.plot(external_date[], label='Données réelles (S)', linestyle='dashed')
    plt.xlabel('Jours')
    plt.ylabel('Nombre de personnes')
    plt.title('Sains (S)')
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(S_prediction, label='E (modèle SEIR)')
    plt.plot(external_date[], label='Données réelles (E)', linestyle='dashed')
    plt.xlabel('Jours')
    plt.ylabel('Nombre de personnes')
    plt.title('Exposés (E)')
    plt.legend()

    plt.tight_layout()
    plt.show()

