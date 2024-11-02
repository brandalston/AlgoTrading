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


def rsi(data, period):
    data['delta'] = data.Close.diff()
    data['gain'] = data.delta.mask(data.delta < 0, 0.0)
    data['loss'] = -data.delta.mask(data.delta > 0, -0.0)
    data['avgGain'] = data.gain.rolling(window=period).mean()
    data['avgLoss'] = data.loss.rolling(window=period).mean()
    data['rs'] = data.avgGain / data.avgLoss
    data['RSI'] = 100 - (100 / (1 + data.rs))
    # trading strat for RSI over trading period
    data['rsiSignal'] = 0
    data.loc[data['RSI'] < 30, 'rsiSignal'] = 1
    data.loc[data['RSI'] > 70, 'rsiSignal'] = -1
    # return of RSI over trading period
    data['rsiReturn'] = data.dailyReturn * data['rsiSignal'].shift(1)
    data['rsiCumulativeReturn'] = (1 + data['rsiReturn']).cumprod()
    return data


def sma(data, period):
    #calc 50-day SMA
    data['SMA'] = data.Close.rolling(window=period).mean()
    data['smaSignal'] = 0
    #trading strat for SMA over trading period
    data['smaSignal'] = np.where(data.Close > data.SMA, 1, 0)
    data['smaReturn'] = data.dailyReturn * data.smaSignal.shift(1)
    data['smaCumulativeReturn'] = (1 + data['smaReturn']).cumprod()
    return data

symbol = 'AAPL'
start_date, end_date = '2012-01-01', '2022-12-31'
data = yf.download(symbol, start=start_date, end=end_date)
data['dailyReturn'] = data.Close.pct_change()

# Calculate the 14-day RSI and 50 day SMA
# data = rsi(data,14)
data = sma(data,50)
data = rsi(data, 14)

# Fetch historical data for SPY

# Calculate daily returns and cumulative returns for SPY
data['spyCumulativeReturn'] = (1 + data['dailyReturn']).cumprod()

# Plot both cumulative returns on the same chart
plt.figure(figsize=(12, 6))
plt.plot(data.index, data['smaCumulativeReturn'], label='SMA Strategy')
plt.plot(data.index, data['rsiCumulativeReturn'], label='RSI Strategy')
plt.plot(data.index, data['spyCumulativeReturn'], label='Reg')
plt.xlabel('Date'); plt.ylabel('Cumulative Returns'); plt.legend()
plt.show()
plt.close()

