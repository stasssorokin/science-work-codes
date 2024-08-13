import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid


Mn = 0.93827
a = 1/137


def calc_sigmaT(channel: str, Q2: float):
    df = pd.read_csv(rf'data\{channel}_sigmaU_data.csv', sep=' ')
    cols = ['Q2', 'W', 'cosTheta', 'sigmaT', 'dsigmaT']
    new_df = pd.DataFrame(columns=cols)
    loc = 0
    for W in sorted(set(df.loc[df['Q2min'] == Q2]['Wmin'])):
        for cosTheta in sorted(set(df.loc[df['Q2min'] == Q2].loc[df['Wmin'] == W]['cosTheta'])):
            sigmaU = df.loc[df['Q2min'] == Q2].loc[df['Wmin'] == W].loc[df['cosTheta'] == cosTheta]['sigmaU'].values[0]
            dsigmaU = df.loc[df['Q2min'] == Q2].loc[df['Wmin'] == W].loc[df['cosTheta'] == cosTheta]['dsigmaU'].values[0]
            if cosTheta > 1:
                cosTheta = round(np.cos(cosTheta / 180 * np.pi), 1)
            nu = (W**2 - Mn**2 + Q2) / (2 * Mn)
            theta = np.arccos(cosTheta)
            tg = np.tan(theta / 2)
            epsilon = (1 + 2 * (1 + nu**2 / Q2) * tg**2)**(-1)
            sigmaT = round(sigmaU / 1.2, 6)
            dsigmaT = round(dsigmaU / 1.2, 6)
            i_data = [Q2, W, cosTheta, sigmaT, dsigmaT]
            new_df.loc[loc] = i_data
            loc +=1
    new_df.to_csv(fr'data\{channel}_sigmaT_calc_data.csv', sep=' ', index=False)


def calc_f2(channel: str):
    cols = ['Q2', 'W', 'F2', 'dF2']
    new_df = pd.DataFrame(columns=cols)
    loc = 0
    df = pd.read_csv(fr'data\\{channel}_sigmaT_calc_data.csv', sep=' ')
    for Q2 in set(df['Q2'].values):
        for W in set(df.loc[df['Q2'] == Q2]['W'].values):
            xx = df.loc[df['Q2'] == Q2].loc[df['W'] == W]['cosTheta'].values
            yy = df.loc[df['Q2'] == Q2].loc[df['W'] == W]['sigmaT'].values
            yy_max = [a + b for a, b in zip(yy, df.loc[df['Q2'] == Q2].loc[df['W'] == W]['dsigmaT'].values)]
            yy_min = [a - b for a, b in zip(yy, df.loc[df['Q2'] == Q2].loc[df['W'] == W]['dsigmaT'].values)]
            sigmaT = 2 * np.pi * abs(trapezoid(y=yy, x=xx)) / 389.5
            dsigmaT = 2 * np.pi * (abs(trapezoid(x=xx, y=yy_max)) - abs(trapezoid(x=xx, y=yy_min))) / 389.5
            K = (W**2 - Mn**2) / (2 * Mn)
            F1 = (K * W * sigmaT) / (4 * np.pi**2 * a)
            dF1 = (K * W * dsigmaT) / (4 * np.pi**2 * a)
            x = Q2 / (W**2 - Mn**2 + Q2)
            F2 = (F1 * 2.4 * x) / (1 + (Q2 * 4 * Mn**2) / (W**2 - Mn**2 + Q2)**2)
            dF2 = (dF1 * 2.4 * x) / (1 + (Q2 * 4 * Mn**2) / (W**2 - Mn**2 + Q2)**2)
            new_df.loc[loc] = [Q2, W, F2, dF2]
            loc += 1
    new_df = new_df.round(6)
    new_df.to_csv(fr'data\{channel}_f2_calc_data.csv', sep=' ', index=False)


def graph(Q2: float):
    df_incl = pd.read_csv(r'data\incl_calc_data.csv', sep=' ')
    df_pi0p = pd.read_csv(r'data\Pi0P_f2_calc_data.csv', sep=' ')
    df_pin = pd.read_csv(r'data\Pin_f2_calc_data.csv', sep=' ')
    fig, axs = plt.subplots(1, 2, figsize=[17, 8], layout='constrained')
    ms = 5
    lfs = 16
    tfs = 16
    cs = 3
    method = 'linear'
    X = df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['W'].values
    axs[0].errorbar(
        x= X,
        y=df_pi0p.loc[df_pi0p['Q2'] == Q2].sort_values('W')['F2'].values[:len(X)],
        yerr=df_pi0p.loc[df_pi0p['Q2'] == Q2].sort_values('W')['dF2'].values[:len(X)],
        fmt='-o', color='red', markersize=ms, capsize=cs,
        label=r'$\gamma_\upsilon p \rightarrow \pi^0 p$')

    axs[0].errorbar(
        x=X,
        y=df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['F2'].values,
        yerr=df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['dF2'].values,
        fmt='-o', color='blue', markersize=ms, capsize=cs, 
        label=r'$\gamma_\upsilon p \rightarrow \pi^+ n$')
    
    xx = df_incl['W'].values
    yy = df_incl['Q2'].values
    points = np.array(list(zip(xx, yy)))
    values = df_incl['F2'].values
    xx_new = X
    yy_new = [Q2 for _ in xx_new]
    xi = np.array(list(zip(xx_new, yy_new)))
    interp_val = griddata(points, values, xi, method=method)

    errs = df_incl['dF2'].values
    interp_errs = griddata(points, errs, xi, method=method)

    axs[0].errorbar(
        x=xx_new, 
        y=interp_val,
        yerr=interp_errs,
        fmt='-o', color='black', markersize=ms, capsize=cs, 
        label='inclusive')

    axs[0].set_xlabel('W, GeV', fontsize=lfs)
    axs[0].set_ylabel(r'F$_2$', fontsize=lfs)
    axs[0].legend(title=rf'Q$^2$ = {Q2} GeV$^2$', title_fontsize=tfs,loc='upper right', fontsize=lfs)
    axs[0].grid(visible=True, linestyle='--')


    y0=df_pi0p.loc[df_pi0p['Q2'] == Q2].sort_values('W')['F2'].values[:len(X)]
    dy0=df_pi0p.loc[df_pi0p['Q2'] == Q2].sort_values('W')['dF2'].values[:len(X)]
    y1=df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['F2'].values
    dy1=df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['dF2'].values
    axs[1].errorbar(
        x= X,
        y=[sum(i) for i in zip(y0, y1)],
        # yerr=[((dy0[i])**2 + (dy1[i])**2)**0.5 for i in range(len(dy0))],
        yerr=[dy0[i] + dy1[i] for i in range(len(dy0))],
        fmt='-o', color='red', markersize=ms, capsize=cs,
        label=r'$\gamma_\upsilon p \rightarrow \pi^0 p \ (\pi^+ n)$')
    
    xx = df_incl['W'].values
    yy = df_incl['Q2'].values
    points = np.array(list(zip(xx, yy)))
    values = df_incl['F2'].values
    xx_new = X
    yy_new = [Q2 for _ in xx_new]
    xi = np.array(list(zip(xx_new, yy_new)))
    interp_val = griddata(points, values, xi, method=method)

    errs = df_incl['dF2'].values
    interp_errs = griddata(points, errs, xi, method=method)

    axs[1].errorbar(
        x=xx_new, 
        y=interp_val,
        yerr=interp_errs,
        fmt='-o', color='black', markersize=ms, capsize=cs, 
        label='inclusive')

    axs[1].set_xlabel('W, GeV', fontsize=lfs)
    axs[1].set_ylabel(r'F$_2$', fontsize=lfs)
    axs[1].legend(title=rf'Q$^2$ = {Q2} GeV$^2$', title_fontsize=tfs, loc='upper right', fontsize=lfs)
    axs[1].grid(visible=True, linestyle='--')

    plt.savefig(rf'graphs\incl_vs_excl_{Q2}.png')
    # plt.show()
    

def graph_frac(Q2: float):
    df_incl = pd.read_csv(r'data\incl_calc_data.csv', sep=' ')
    df_pi0p = pd.read_csv(r'data\Pi0P_f2_calc_data.csv', sep=' ')
    df_pin = pd.read_csv(r'data\Pin_f2_calc_data.csv', sep=' ')
    fig, ax = plt.subplots(figsize=[17, 8], layout='constrained')
    ms = 5
    lfs = 16
    tfs = 16
    cs = 3
    method = 'linear'
    X = df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['W'].values

    y0=df_pi0p.loc[df_pi0p['Q2'] == Q2].sort_values('W')['F2'].values[:len(X)]
    dy0=df_pi0p.loc[df_pi0p['Q2'] == Q2].sort_values('W')['dF2'].values[:len(X)]
    y1=df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['F2'].values
    dy1=df_pin.loc[df_pin['Q2'] == Q2].sort_values('W')['dF2'].values
    y2=[sum(i) for i in zip(y0, y1)]
    dy2=[((dy0[i])**2 + (dy1[i])**2)**0.5 for i in range(len(dy0))]
    
    xx = df_incl['W'].values
    yy = df_incl['Q2'].values
    points = np.array(list(zip(xx, yy)))
    values = df_incl['F2'].values
    xx_new = X
    yy_new = [Q2 for _ in xx_new]
    xi = np.array(list(zip(xx_new, yy_new)))
    interp_val = griddata(points, values, xi, method=method)

    errs = df_incl['dF2'].values
    interp_errs = griddata(points, errs, xi, method=method)

    y = [a / b for a, b in zip(interp_val, y2)]
    yerr=[(a * db + b * da)/b**2 for a, da, b, db in zip(interp_val, interp_errs, y2, dy2)]
    ax.errorbar(
        x=xx_new, 
        y=y,
        # yerr=[(a**2 + b**2)**0.5 for a, b in zip(interp_errs, dy2)],
        yerr=yerr,
        fmt='-o', color='black', markersize=ms, capsize=cs, 
        label='inclusive')

    ax.set_xlabel('W, GeV', fontsize=lfs)
    ax.set_ylabel(r'F$_2$(incl) / F$_2$(sum)', fontsize=lfs)
    ax.legend(title=rf'Q$^2$ = {Q2} GeV$^2$', title_fontsize=tfs, loc='upper right', fontsize=lfs)
    ax.grid(visible=True, linestyle='--')

    plt.savefig(rf'graphs\frac_incl_excl_{Q2}.png')
    # plt.show()


Q2 = 0.35

calc_sigmaT(channel='Pi0P', Q2=Q2)
calc_sigmaT(channel='Pin', Q2=Q2)

calc_f2('Pi0P')
calc_f2('Pin')

graph(Q2)

graph_frac(Q2)
