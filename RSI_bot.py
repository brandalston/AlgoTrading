import numpy as np, pandas as pd, yfinance as yf, matplotlib.pyplot as plt

symbol = 'QQQ'
start_date, end_date =  '2015-01-01', '2022-12-31'
data = yf.download(symbol, start=start_date, end=end_date)

def rsi(data, period):
    delta = data.diff().dropna()
    gain = delta.where(delta > 0, 0)
    loss = -delta.where(delta < 0, 0)
    avg_gain = gain.rolling(window=period).mean()
    avg_loss = loss.rolling(window=period).mean()
    rs = avg_gain / avg_loss
    return 100 - (100 / (1 + rs))

# Calculate the 14-day RSI
data['RSI'] = rsi(data['Close'], 14)

#trading strat for RSI over trading period
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
plt.plot(data.index, data['Cumulative_Return'], label='RSI Strategy')
plt.plot(spy_data.index, spy_data['Cumulative_Return'], label='SPY')
plt.xlabel('Date'); plt.ylabel('Cumulative Returns'); plt.legend()
plt.show(); plt.close()
