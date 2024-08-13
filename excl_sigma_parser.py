import requests
from bs4 import BeautifulSoup
import re
import pandas as pd


url = 'https://clas.sinp.msu.ru/cgi-bin/almaz/form1.py'

def get_data(channel: str):

    payload = {
        'Particle': channel,
        'Q2min': '0.25',
        'Q2max': '5',
        'Wmin': '1.11',
        'Wmax': '2'
    }

    response = requests.get(url, params=payload)

    if response.ok:
        print(channel)
        cols = ['Q2min', 'Q2max', 'Wmin', 'Wmax', 'E', 'cosTheta', 'sigmaT', 'dsigmaT', 'sigmaL', 'dsigmaL', 'sigmaTT', 'dsigmaTT', 'sigmaLT', 'dsigmaLT']
        df = pd.DataFrame(columns = cols)
        loc = 0
        url_loc = 1
        pattern = r'Q2min=(?P<Q2min>\d?\.?\d+)\&Q2max=(?P<Q2max>\d?\.?\d+)\&Wmin=(?P<Wmin>\d?\.?\d+)\&Wmax=(?P<Wmax>\d?\.?\d+)\&E=(?P<E>\d?\.?\d+)'
        soup = BeautifulSoup(response.text, 'lxml')
        all_a = soup.find_all('a')
        url_num = len(all_a)
        urls = []
        for i_a in all_a:
            urls.append('https://clas.sinp.msu.ru/cgi-bin/almaz/' + i_a.get('href'))

        for i_url in urls:
            i_response = requests.get(i_url)
            if i_response.ok:
                m = re.search(pattern, i_url)
                tds = [m['Q2min'], m['Q2max'], m['Wmin'], m['Wmax'], m['E']]
                i_soup = BeautifulSoup(i_response.text, 'lxml')
                all_tr = i_soup.find_all('tr')
                for i in range(1, len(all_tr)):
                    for i_td in all_tr[i].find_all('td'):
                        tds.append(i_td.text)
                    df.loc[loc] = tds
                    loc += 1
                    tds = [m['Q2min'], m['Q2max'], m['Wmin'], m['Wmax'], m['E']]
            else: print(f'Error for {i_url}')
            print(f'Progress: {url_loc}/{url_num}')
            url_loc += 1
        df.to_csv(f'diploma\\data\\{channel.lower()}_sigma_data.csv', sep=' ', index=False)

    else: print(f'Error for {url}')

get_data('Pi0P')
get_data('Pin')