import re
import pandas as pd
import math


def get_pi0p_data():
    with open('diploma\\data\\pi0p_sigma_clasdb_data.txt') as file:
        text = file.read()
        lines = text.split('\n')
        pattern_Q2 = r'Q\^2 range: (?P<Q2min>\d+\.?\d*) - (?P<Q2max>\d+\.?\d*)'
        pattern_sigma = r'(?P<W>\d+\.?\d*)\s*(?P<epsilon>\d+\.?\d*)\s*(?P<cosTheta>-?\d+\.?\d*)\s*(?P<phi>\d+\.?\d*)\s*(?P<sigma>\d+\.?\d*)\s*(?P<dsigma_stat>\d+\.?\d*)\s*(?P<dsigma_sys>\d+\.?\d*)'
        df = pd.DataFrame(columns = ['Q2min', 'Q2max', 'beam_energy', 'W', 'epsilon', 'cosTheta', 'phi', 'sigma', 'dsigma_stat', 'dsigma_sys', 'dsigma'])
        loc = 0
        for line in lines:
            m = re.search(pattern_Q2, line)
            if m:
                Q2min = float(m['Q2min'])
                Q2max = float(m['Q2max'])
                continue
            m = re.search(pattern_sigma, line)
            if m:
                if Q2min < 1.3:
                    beam_energy = 1.645
                else:
                    beam_energy = 2.445
                dsigma = float(m['dsigma_stat']) + float(m['dsigma_sys'])
                i_data = [Q2min, Q2max, beam_energy] + list(m.groups()) + [dsigma]
                df.loc[loc] = i_data
                loc += 1
        else:
            df = df.round(6)
            df = df[['Q2min', 'Q2max', 'W', 'epsilon', 'beam_energy', 'cosTheta', 'phi', 'sigma', 'dsigma_stat', 'dsigma_sys', 'dsigma']]
            df.to_csv('diploma\\data\\pi0p_sigma_clasdb_data.csv', sep=' ', index=False)

def get_pin_data():
    with open('diploma\\data\\pin_sigma_clasdb_data.txt') as file:
        text = file.read()
        lines = text.split('\n')
        pattern_Q2 = r'Q\^2 range: (?P<Q2min>\d+\.?\d*) - (?P<Q2max>\d+\.?\d*)'
        pattern_W = r'W range: (?P<Wmin>\d+\.?\d*) - (?P<Wmax>\d+\.?\d*)'
        pattern_sigma1 = r'^(?P<beam_energy>\d+\.?\d*)\s*(?P<theta>\d+\.?\d*)\s*(?P<phi>\d+\.?\d*)\s*(?P<sigma>\d+\.?\d*)\s*(?P<dsigma>\d+\.?\d*)$'
        pattern_sigma2 = r'^(?P<W>\d+\.?\d*)\s*(?P<Q2>\d+\.?\d*)\s*(?P<beam_energy>\d+\.?\d*)\s*(?P<epsilon>\d+\.?\d*)\s*(?P<cosTheta>-?\d+\.?\d*)\s*(?P<phi>\d+\.?\d*)\s*(?P<sigma>\d+\.?\d*)\s*(?P<dsigma_stat>\d+\.?\d*)\s*(?P<dsigma_sys>\d+\.?\d*)$'
        cols = ['Q2min', 'Q2max', 'Wmin', 'Wmax', 'beam_energy', 'cosTheta', 
                'phi', 'sigma', 'dsigma_stat', 'dsigma_sys', 'dsigma']
        df = pd.DataFrame(columns=cols)
        loc = 0
        for line in lines:
            m = re.search(pattern_Q2, line)
            if m:
                Q2min = m['Q2min']
                Q2max = m['Q2max']
                continue
            m = re.search(pattern_W, line)
            if m:
                Wmin = m['Wmin']
                Wmax = m['Wmax']
                continue
            m = re.search(pattern_sigma1, line)
            if m:
                theta = float(m['theta']) * math.pi / 180
                cosTheta = (math.cos(theta))
                dsigma_stat = m['dsigma']
                dsigma_sys = 0
                dsigma = m['dsigma']
                i_data = [Q2min, Q2max, Wmin, Wmax, m['beam_energy'], cosTheta, 
                          m['phi'], m['sigma'], dsigma_stat, dsigma_sys, dsigma]
                df.loc[loc] = i_data
                loc += 1
                continue
            m = re.search(pattern_sigma2, line)
            if m:
                Q2min = m['Q2']
                Q2max = m['Q2']
                Wmin = m['W']
                Wmax = m['W']
                cosTheta = m['cosTheta']
                dsigma_stat = m['dsigma_stat']
                dsigma_sys = m['dsigma_sys']
                dsigma = (float(dsigma_stat) + float(dsigma_sys))
                i_data = [Q2min, Q2max, Wmin, Wmax, m['beam_energy'], cosTheta, 
                          m['phi'], m['sigma'], dsigma_stat, dsigma_sys, dsigma]
                df.loc[loc] = i_data
                loc += 1
        df = df.astype(float)
        df = df.round(6)
        df.to_csv('diploma\\data\\pin_sigma_clasdb_data.csv', sep=' ', index=False)

# get_pi0p_data()
# get_pin_data()
