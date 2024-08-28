import numpy as np
import pandas as pd
from tabulate import tabulate
time_call_prices = []
time_put_prices = []

S0,K,T,r,sigma,M=100,105,5,0.05,0.4,20
times2 = [0, 0.5, 1, 1.5, 3, 4.5]			
call_prices = []
put_prices = []

def n_step_pricing(fin_pricing, fin_pricing2, M, p, q, rate):	
	if(M==0 or M==2 or M==4 or M==6 or M==12 or M==18):
		time_call_prices.append(np.round(fin_pricing,3))
		time_put_prices.append(np.round(fin_pricing2,3))
	if(M==0):
		return fin_pricing[0], fin_pricing2[0]
	fin_pricing_step = []	
	fin_pricing_step2 = []
	for i in range(len(fin_pricing)-1):
		fin_pricing_step.append((p*fin_pricing[i] + q*fin_pricing[i+1])/rate)
		fin_pricing_step2.append((p*fin_pricing2[i] + q*fin_pricing2[i+1])/rate)
	return n_step_pricing(fin_pricing_step, fin_pricing_step2, M-1, p, q, rate)


u = np.exp(sigma*np.sqrt(T/M) + (r-sigma*sigma/2)*T/M)
d = np.exp(-sigma*np.sqrt(T/M) + (r-sigma*sigma/2)*T/M)
rate = np.exp(r*T/M)
p = (rate-d)/(u-d)
q = 1-p

if(p<0 or p>1):
	exit("Arbitrage exists")

stock_prices = []

for i in range(M+1):
	stock_prices.append(S0*np.math.pow(u,M-i)*np.math.pow(d,i))

fin_pricing = []
fin_pricing2 = []

for i in range(M+1):
	fin_pricing.append(max(stock_prices[i]-K, 0))
	fin_pricing2.append(max(K-stock_prices[i], 0))

call_price, put_price= n_step_pricing(fin_pricing, fin_pricing2, M, p, q, rate)
call_prices.append(call_price)
put_prices.append(put_price)
df1 = pd.DataFrame({"Time":times2, "Call Price":time_call_prices[::-1]})
df2 = pd.DataFrame({"Time":times2, "Put Price":time_put_prices[::-1]})
print(tabulate(df1, headers='keys', tablefmt='grid', showindex=False))
print(tabulate(df2, headers='keys', tablefmt='grid', showindex=False))