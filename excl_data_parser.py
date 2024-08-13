from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import re
import pandas as pd
import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import os


TIMEOUT = 30
Mn = 0.93827
ALPHA = 1 / 137

def get_source_links(channel: str):
    URL = 'https://clas.sinp.msu.ru/~almaz/EVALUATED_EXCLUSIVE_STRUCTURE_FUNCTIONS/FolderY/index.html'
    log_file = open(f'{channel}_source.log', '+a')
    log_file.truncate(0)
    data_file = open(f'{channel}_source_links.txt', '+a')
    data_file.truncate(0)
    chrome_options = Options()
    chrome_options.add_argument('--headless')
    driver = webdriver.Chrome(options=chrome_options)
    print(channel)
    print(channel, file=log_file)
    driver.get(URL)
    Select(driver.find_element(By.NAME, 'Particle')).select_by_value(channel)
    driver.find_element(By.NAME, 'Q2min').send_keys('0.2')
    driver.find_element(By.NAME, 'Q2max').send_keys('5')
    driver.find_element(By.NAME, 'Wmin').send_keys('1.11')
    driver.find_element(By.NAME, 'Wmax').send_keys('2.5')
    driver.find_element(By.CLASS_NAME, 'button').click()
    url_list = [i_url.get_attribute('href') for i_url in driver.find_elements(By.TAG_NAME, 'a')]
    n_total = len(url_list)
    url_loc = 1
    for i_url in url_list:
        print(f'{url_loc}/{n_total}')
        print(f'{url_loc}/{n_total}', file=log_file)
        url_loc += 1
        print(i_url)
        print(i_url, file=log_file)
        try:
            driver.get(i_url)
            WebDriverWait(driver, TIMEOUT).until(EC.visibility_of_element_located((By.TAG_NAME, 'a')))
            source_link = driver.find_element(By.TAG_NAME, 'a').get_attribute('href')
            print(source_link)
            print(source_link, file=log_file)
            data_file.write(source_link + '\n')
        except:
            print('ERROR:', i_url)
            print('ERROR:', i_url, file=log_file)
    driver.quit()
    log_file.close()
    data_file.close()


def get_sigmaU_data(channel: str):
    df = pd.DataFrame(columns=['Q2min', 'Q2max', 'Q2', 'Wmin', 'Wmax', 'W', 'cosTheta', 'epsilon', 'sigmaU', 'dsigmaU'])
    loc = 0
    file = open(f'{channel}_source_links.txt', 'r')
    log_file = open(f'{channel}_sigmaU_data1.log', '+a')
    log_file.truncate(0)
    url_list = file.read().splitlines()
    if channel == 'Pin':
        url_list = url_list[:258]
    url_loc = 1
    url_total = len(url_list)
    Q2_pattern1 = r'Q2 : (?P<Q2min>\d+\.?\d*) . (?P<Q2max>\d+\.?\d*)'
    Q2_pattern2 = r'Q2 : (?P<Q2>\d+\.?\d*)'
    W_pattern1 = r'W : (?P<Wmin>\d+\.?\d*) . (?P<Wmax>\d+\.?\d*)'
    W_pattern2 = r'W : (?P<W>\d+\.?\d*)'
    epsilon_pattern = r'ε : (?P<epsilon>\d+\.?\d*)'
    beam_energy_pattern = r'Ebeam : (?P<beam_energy>\d+\.?\d*)'
    chrome_options = Options()
    chrome_options.add_argument('--headless')
    driver = webdriver.Chrome(options=chrome_options)
    for i_url in url_list:
        print(f'{url_loc}/{url_total}')
        print(f'{url_loc}/{url_total}', file=log_file)
        url_loc +=1
        print(i_url)
        print(i_url, file=log_file)
        driver.get(i_url)
        WebDriverWait(driver, TIMEOUT).until(EC.visibility_of_element_located((By.CLASS_NAME, 'props')))
        props = driver.find_element(By.CLASS_NAME, 'props').text.splitlines()
        m = re.search(pattern=Q2_pattern1, string=props[0])
        if m:
            Q2min = float(m['Q2min'])
            Q2max = float(m['Q2max'])
        else:
            m = re.search(pattern=Q2_pattern2, string=props[0])
            Q2min = float(m['Q2'])
            Q2max = Q2min
        if Q2min == 0.35:
            Q2 = 0.4
        elif Q2min == 0.45:
            Q2 = 0.5
        elif (Q2min == 0.55) or (Q2min == 0.6):
            Q2 = 0.625
        else:
            Q2 = round((Q2min + Q2max) / 2, 2)
        m = re.search(pattern=W_pattern1, string=props[1])
        if m:
            Wmin = float(m['Wmin'])
            Wmax = float(m['Wmax'])
        else:
            m = re.search(pattern=W_pattern2, string=props[1])
            Wmin = float(m['W'])
            Wmax = Wmin
        W = round((Wmin + Wmax) / 2, 2)
        m = re.search(pattern=epsilon_pattern, string=props[2])
        if m:
            epsilon = float(m['epsilon'])
        else:
            m = re.search(pattern=beam_energy_pattern, string=props[2])
            beam_energy = float(m['beam_energy'])
            nu = (W**2 - Mn**2 + Q2) / (2 * Mn)
            sinTheta2e = (Q2 / (4*beam_energy*(beam_energy - nu)))**0.5
            epsilon = round((1 + 2*(1 + nu**2/Q2)*(np.tan(np.arcsin(sinTheta2e)))**2)**(-1), 3)
        WebDriverWait(driver, TIMEOUT).until(EC.visibility_of_element_located((By.NAME, 'fit1')))
        driver.find_element(By.NAME, 'fit1').click()
        WebDriverWait(driver, TIMEOUT).until(EC.visibility_of_element_located((By.NAME, 'fixval')))
        option_list = Select(driver.find_element(By.NAME, 'fixval')).options
        for option in option_list: # fixed cosTheta
            option.click()
            cosTheta = float(option.text)
            if cosTheta > 1:
                cosTheta = round(np.cos(cosTheta / 180 * np.pi), 2)
            WebDriverWait(driver, TIMEOUT).until(EC.visibility_of_element_located((By.TAG_NAME, 'iframe')))
            driver.switch_to.frame("coeff1")
            WebDriverWait(driver, TIMEOUT).until(EC.visibility_of_all_elements_located((By.TAG_NAME, 'tr')))
            tr_list = driver.find_elements(By.TAG_NAME, 'tr')
            for tr in tr_list[2:3]: # fixed rows
                WebDriverWait(driver, TIMEOUT).until(EC.visibility_of_all_elements_located((By.TAG_NAME, 'td')))
                td_list = tr.find_elements(By.TAG_NAME, 'td')
                sigmaU = float(td_list[0].text)
                dsigmaU = float(td_list[1].text)
                i_data = [Q2min, Q2max, Q2, Wmin, Wmax, W, cosTheta, epsilon, sigmaU, dsigmaU]
                print(i_data)
                print(i_data, file=log_file)
                df.loc[loc] = i_data
                loc +=1
            driver.switch_to.default_content()
    df.to_csv(fr'data\{channel}_sigmaU_data1.csv', sep=' ', index=False)
    driver.quit()
    file.close()
    log_file.close()


def graph_sigmaU(channel: str):
    df = pd.read_csv(fr'data\{channel}_sigmaU_data1.csv', sep=' ')
    # Q2min_set = sorted(set(df['Q2min']))
    # Q2max_set = sorted(set(df['Q2max']))
    Q2_set = sorted(set(df['Q2']))
    if channel == 'Pi0P':
        channel_label = r'$\gamma_\upsilon p \rightarrow \pi^0 p$'
    else:
        channel_label = r'$\gamma_\upsilon p \rightarrow \pi^+ n$'
    for Q2 in Q2_set: # Q2min_set
        # Wmin_set = sorted(set(df.loc[df['Q2min'] == Q2]['Wmin']))
        # Wmax_set = sorted(set(df.loc[df['Q2min'] == Q2]['Wmax']))
        W_set = sorted(set(df.loc[df['Q2'] == Q2]['W']))
        for W in W_set: # Wmin_set
            epsilon_set = sorted(set(df.loc[df['Q2'] == Q2].loc[df['W'] == W]['epsilon']))
            for epsilon in epsilon_set:
                print(f'{channel:<5}{Q2:<5}{W:<5}{epsilon:<5}')
                # x_data = df.loc[df['Q2min'] == Q2].loc[df['Wmin'] == W].sort_values('cosTheta')['cosTheta'].values
                # y_data = df.loc[df['Q2min'] == Q2].loc[df['Wmin'] == W].sort_values('cosTheta')['sigmaU'].values
                # yerr_data = df.loc[df['Q2min'] == Q2].loc[df['Wmin'] == W].sort_values('cosTheta')['dsigmaU'].values
                x_data = df.loc[df['Q2'] == Q2].loc[df['W'] == W].loc[df['epsilon'] == epsilon].sort_values('cosTheta')['cosTheta'].values
                y_data = df.loc[df['Q2'] == Q2].loc[df['W'] == W].loc[df['epsilon'] == epsilon].sort_values('cosTheta')['sigmaU'].values
                yerr_data = df.loc[df['Q2'] == Q2].loc[df['W'] == W].loc[df['epsilon'] == epsilon].sort_values('cosTheta')['dsigmaU'].values
                
                plt.figure(figsize=(11, 9))
                plt.errorbar(x=x_data, y=y_data, yerr=yerr_data, fmt='ok', markersize=5, capsize=3, 
                            label=channel_label + '\n' + fr'Q$^2$ = {Q2} GeV$^2$' + '\n' + f'W = {W} GeV' + '\n' + fr'$\varepsilon$ = {epsilon}')
                plt.xlabel(r'$\cos(\theta)$', fontsize=18)
                plt.ylabel(r'$\frac{d\sigma_U}{d\cos(\theta)}$, mcb', fontsize=16)
                plt.legend(loc='upper right', fontsize=16)
                plt.grid(visible=True, linestyle='--')
            
                path = os.path.dirname(__file__) + fr'\graphs1\{channel}\Q2\{Q2}'
                if not os.path.isdir(path):
                    os.makedirs(path)
                plt.savefig(path + fr'\{channel}_sigmaU_Q2_{Q2}_W_{W}_eps_{epsilon}.png')
                
                path = os.path.dirname(__file__) + fr'\graphs1\{channel}\W\{W}'
                if not os.path.isdir(path):
                    os.makedirs(path)
                plt.savefig(path + fr'\{channel}_sigmaU_W_{W}_Q2_{Q2}_eps_{epsilon}.png')
                plt.close()


def calc_f2_data(channel: str):
    print(channel)
    new_df = pd.DataFrame(columns=['Q2', 'W', 'epsilon', 'F2', 'dF2'])
    loc = 0
    RLT = 0.8
    df = pd.read_csv(fr'data\{channel}_sigmaU_data1.csv', sep=' ')
    df['sigmaT'] = round(df['sigmaU'] / (1 + RLT*df['epsilon']), 4)
    df['dsigmaT'] = round(df['dsigmaU'] / (1 + RLT*df['epsilon']), 4)
    # Q2_set = sorted(set(df['Q2']))
    Q2_set = (0.4, 0.5, 0.625)
    for Q2 in Q2_set: # Q2min_set[:1]
        W_set = sorted(set(df.loc[df['Q2'] == Q2]['W']))
        for W in W_set: # Wmin_set[:1]
            epsilon_set = sorted(set(df.loc[df['Q2'] == Q2].loc[df['W'] == W]['epsilon']))
            for epsilon in epsilon_set:
                x_data = df.loc[df['Q2'] == Q2].loc[df['W'] == W].loc[df['epsilon'] == epsilon].sort_values('cosTheta')['cosTheta'].values
                y_data = df.loc[df['Q2'] == Q2].loc[df['W'] == W].loc[df['epsilon'] == epsilon].sort_values('cosTheta')['sigmaT'].values
                y_data_max = [round(a + b, 4) for a, b in zip(y_data, 
                    df.loc[df['Q2'] == Q2].loc[df['W'] == W].loc[df['epsilon'] == epsilon].sort_values('cosTheta')['dsigmaT'].values)]
                y_data_min = [round(a - b, 4) for a, b in zip(y_data, 
                    df.loc[df['Q2'] == Q2].loc[df['W'] == W].loc[df['epsilon'] == epsilon].sort_values('cosTheta')['dsigmaT'].values)]
                sigmaT = round(2 * np.pi * abs(trapezoid(y=y_data, x=x_data)) / 389.5, 4)
                dsigmaT = round(2 * np.pi * abs(trapezoid(y=y_data_max, x=x_data) - trapezoid(y=y_data_min, x=x_data)) / 389.5, 4)
                K = (W**2 - Mn**2) / (2 * Mn)
                nu = (W**2 - Mn**2 + Q2) / (2 * Mn)
                x = Q2 / (2 * Mn * nu)
                F2 = round((K * W * 2 * x * (1 + RLT) * sigmaT) / (4 * (np.pi)**2 * ALPHA * (1 + Q2 / nu**2)), 6)
                dF2 = round((K * W * 2 * x * (1 + RLT) * dsigmaT) / (4 * (np.pi)**2 * ALPHA * (1 + Q2 / nu**2)), 6)
                i_data = [Q2, W, epsilon, F2, dF2]
                new_df.loc[loc] = i_data
                loc += 1
                print(i_data)
    new_df.to_csv(fr'data\{channel}_F2_calc_data1.csv', sep=' ', index=False)


def graph_f2(Q2: float):
    df_incl = pd.read_csv(r'data\incl_calc_data.csv', sep=' ')
    df_pi0p = pd.read_csv(r'data\Pi0P_F2_calc_data1.csv', sep=' ')
    df_pin = pd.read_csv(r'data\Pin_F2_calc_data1.csv', sep=' ')
    fig, axs = plt.subplots(1, 2, figsize=[17, 8], layout='constrained')
    ms = 5
    lfs = 16
    tfs = 16
    cs = 3
    method = 'linear'
    x_data_pi0p = df_pi0p.loc[df_pi0p['Q2'] == Q2].sort_values('W')['W'].values
    x_data_pin = df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['W'].values
    for W in x_data_pi0p:
        eps_pi0p = df_pi0p.loc[df_pi0p['Q2'] == Q2].loc[df_pi0p['W'] == W]['epsilon'].values
        if len(eps_pi0p) > 1:
            df_pi0p = df_pi0p.drop(df_pi0p.loc[df_pi0p['Q2'] == Q2].loc[df_pi0p['W'] == W].loc[df_pi0p['epsilon'] == eps_pi0p[1]].index)
    for W in x_data_pin:
        eps_pin = df_pi0p.loc[df_pi0p['Q2'] == Q2].loc[df_pi0p['W'] == W]['epsilon'].values
        if len(eps_pin) > 1:
            df_pin = df_pin.drop(df_pin.loc[df_pin['Q2'] == Q2].loc[df_pin['W'] == W].loc[df_pin['epsilon'] == eps_pin[1]].index)
    x_data_pin = df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['W'].values
    y_data_pin = df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['F2'].values
    dy_data_pin = df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['dF2'].values
    x_data_pi0p = df_pi0p.loc[df_pi0p['Q2'] == Q2].loc[df_pi0p['W'] < x_data_pin[-1]].sort_values('W')['W'].values
    y_data_pi0p = df_pi0p.loc[df_pi0p['Q2'] == Q2].loc[df_pi0p['W'] < x_data_pin[-1]].sort_values('W')['F2'].values
    dy_data_pi0p = df_pi0p.loc[df_pi0p['Q2'] == Q2].loc[df_pi0p['W'] < x_data_pin[-1]].sort_values('W')['dF2'].values
    axs[1].errorbar(
        x=x_data_pi0p,
        y=y_data_pi0p,
        yerr = dy_data_pi0p,
        fmt='-o', color='red', markersize=ms, capsize=cs,
        label=r'$\gamma_\upsilon p \rightarrow \pi^0 p$'
    )
    axs[1].errorbar(
        x=x_data_pin,
        y = y_data_pin,
        yerr = dy_data_pin,
        fmt='-o', color='blue', markersize=ms, capsize=cs, 
        label=r'$\gamma_\upsilon p \rightarrow \pi^+ n$'
    )
    x_data_incl = df_incl['W'].values
    y_data_incl = df_incl['Q2'].values
    points = np.array(list(zip(x_data_incl, y_data_incl)))
    values = df_incl['F2'].values
    x_new = x_data_pi0p
    y_new = [Q2 for _ in x_new]
    xi = np.array(list(zip(x_new, y_new)))
    interp_val = griddata(points, values, xi, method=method)
    errs = df_incl['dF2'].values
    interp_errs = griddata(points, errs, xi, method=method)
    axs[1].errorbar(
        x=x_new, 
        y=interp_val,
        yerr=interp_errs,
        fmt='-o', color='black', markersize=ms, capsize=cs, 
        label='inclusive'
    )
    axs[1].set_xlabel('W, GeV', fontsize=lfs)
    axs[1].set_ylabel(r'F$_2$', fontsize=lfs)
    axs[1].legend(title=rf'Q$^2$ = {Q2} GeV$^2$', title_fontsize=tfs,loc='upper right', fontsize=lfs)
    axs[1].grid(visible=True, linestyle='--')
    df_pin['W'] = df_pin['W'] - 0.01
    x_data = df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['W'].values
    y_data_pi0p = df_pi0p.loc[df_pi0p['Q2'] == Q2].loc[df_pi0p['W'] <= x_data[-1]].sort_values('W')['F2'].values
    y_data = [a + b for a, b in zip(y_data_pi0p, y_data_pin)]
    dy_data_pi0p = df_pi0p.loc[df_pi0p['Q2'] == Q2].loc[df_pi0p['W'] <= x_data[-1]].sort_values('W')['dF2'].values
    dy_data = [a + b for a, b in zip(dy_data_pi0p, dy_data_pin)]
    axs[0].errorbar(
        x=x_data,
        y=y_data,
        yerr=dy_data,
        fmt='-o', color='red', markersize=ms, capsize=cs,
        label=r'$\gamma_\upsilon p \rightarrow \pi^0 p \ (\pi^+ n)$'
    )
    x_new = x_data
    y_new = [Q2 for _ in x_new]
    xi = np.array(list(zip(x_new, y_new)))
    interp_val = griddata(points, values, xi, method=method)
    errs = df_incl['dF2'].values
    interp_errs = griddata(points, errs, xi, method=method)
    axs[0].errorbar(
        x=x_new, 
        y=interp_val,
        yerr=interp_errs,
        fmt='-o', color='black', markersize=ms, capsize=cs, 
        label='inclusive'
    )
    axs[0].set_xlabel('W, GeV', fontsize=lfs)
    axs[0].set_ylabel(r'F$_2$', fontsize=lfs)
    axs[0].legend(title=rf'Q$^2$ = {Q2} GeV$^2$', title_fontsize=tfs,loc='upper right', fontsize=lfs)
    axs[0].grid(visible=True, linestyle='--')
    path = os.path.dirname(__file__) + fr'\graphs1\comparsion1'
    if not os.path.isdir(path):
        os.makedirs(path)
    plt.savefig(path + fr'\F2_Q2_{Q2}.png')
    plt.close()

# get_source_links(channel='Pi0P')
# get_source_links(channel='Pin')

# get_sigmaU_data('Pi0P')
# get_sigmaU_data('Pin')

# graph_sigmaU(channel='Pi0P')
# graph_sigmaU(channel='Pin')

# calc_f2_data(channel='Pi0P')
# calc_f2_data(channel='Pin')

# Q2_set = [0.4, 0.5, 0.625]
# # graph_f2(0.5)
# for Q2 in Q2_set:
#     graph_f2(Q2)
