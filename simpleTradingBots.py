import numpy as np, pandas as pd, yfinance as yf, matplotlib.pyplot as plt

import requests
import urllib.request
from html.parser import HTMLParser

class GridParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.grid = {}
        self.max_x, self.max_y = 0, 0
        self.in_td = False
        self.current_data = []

    def handle_starttag(self,tag,attrs):
        if tag=='td':
            self.in_td=True

    def handle_endtag(self,tag):
        if tag=='td':
            self.in_td=False
        elif tag=='tr':
            if len(self.current_data)==3:
                x, char, y = self.current_data
                try:
                    x, y = int(x), int(y)
                    self.grid[(x,y)]=char
                    self.max_x = max(self.max_x, x)
                    self.max_y = max(self.max_y, y)
                except ValueError:
                    pass
            self.current_data=[]

    def handle_data(self,data):
        if self.in_td:
            self.current_data.append(data.strip())

def print_grid_from_url(url):
    try:
        print(f"Fecthing URL: {url}")
        with urllib.request.urlopen(url) as response:
            html = response.read().decode('utf-8')
        print(f"Fectched {len(html)} characters of HTML")
        parser = GridParser()
        parser.feed(html)
        print(f"Grid size: {parser.max_x+1} x {parser.max_y+1}")
        for y in range(parser.max_y+1):
            row=''.join(parser.grid.get((x,y),'')for x in range(parser.max_x+1))
            print(row)
    except:
        print('Error parsing url')
url="https://docs.google.com/document/d/e/2PACX-1vQGUck9HIFCyezsrBSnmENk5ieJuYwpt7YHYEzeNJkIb9OSDdx-ov2nRNReKQyey-cwJOoEKUhLmN9z/pub"
url1="https://docs.google.com/document/d/e/2PACX-1vSHesOf9hv2sPOntssYrEdubmMQm8lwjfwv6NPjjmIRYs_FOYXtqrYgjh85jBUebK9swPXh_a5TJ5Kl/pub"

symbol = 'QQQ'
start_date, end_date = '2015-01-01', '2022-12-31'
data = yf.download(symbol, start=start_date, end=end_date)


def rsi(data, period):
    delta = data.diff().dropna()
    gain = delta.where(delta > 0, 0)
    loss = -delta.where(delta < 0, 0)
    avg_gain = gain.rolling(window=period).mean()
    avg_loss = loss.rolling(window=period).mean()
    rs = avg_gain / avg_loss
    return 100 - (100 / (1 + rs))

def sma(data):
    #calc 50-day SMA
    data['SMA_50'] = data['Close'].rolling(window=50).mean()
    
    #trading strat for SMA over trading period
    data['SMA_Signal'] = np.where(data['Close'] > data['SMA_50'], 1, 0)
    data['SMA_Daily_Return'] = data['Close'].pct_change()
    data['SMA_Return'] = data['Daily_Return'] * data['SMA_Signal'].shift(1)
    data['SMA_Cumulative_Return'] = (1 + data['SMA_Return']).cumprod()
    return data

# Calculate the 14-day RSI
data['RSI'] = rsi(data['Close'], 14)
data = sma(data)
# trading strat for RSI over trading period
data['Signal'] = 0
data.loc[data['RSI'] < 30, 'Signal'] = 1
data.loc[data['RSI'] > 70, 'Signal'] = -1

data['Daily_Return'] = data['Close'].pct_change()
data['Strategy_Return'] = data['Daily_Return'] * data['Signal'].shift(1)
data['Cumulative_Return'] = (1 + data['Strategy_Return']).cumprod()

# Fetch historical data for SPY
spy_data = yf.download('SPY', start=start_date, end=end_date)

# Calculate daily returns and cumulative returns for SPY
spy_data['Daily_Return'] = spy_data['Close'].pct_change()
spy_data['Cumulative_Return'] = (1 + spy_data['Daily_Return']).cumprod()

# Plot both cumulative returns on the same chart
plt.figure(figsize=(12, 6))
plt.plot(data.index, data['Cumulative_Return'], label='SMA Strategy')
plt.plot(spy_data.index, spy_data['Cumulative_Return'], label='SPY')
plt.xlabel('Date'); plt.ylabel('Cumulative Returns'); plt.legend()
plt.close()

