import numpy as np
import pandas as pd
def n_step_pricing(fin_pricing, M, p, q, rate):
	if(M==0):
		return fin_pricing[0]
	
	fin_pricing_step = []		
	
	for i in range(len(fin_pricing)-1):
		fin_pricing_step.append((p*fin_pricing[i] + q*fin_pricing[i+1])/rate)
	
	return n_step_pricing(fin_pricing_step, M-1, p, q, rate)	

S0,K,T,r,sigma=100,105,5,0.05,0.4
step = [1,5,10,20,50,100,200,400]
call_prices = []
put_prices = []
for M in step:
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
	
	for i in range(M+1):
		fin_pricing.append(max(stock_prices[i]-K, 0))
	
	call_price = n_step_pricing(fin_pricing, M, p, q, rate)
	put_price = call_price + K*np.exp(-r*T) - S0
	call_prices.append(call_price)
	put_prices.append(put_price)

df = pd.DataFrame({"M":step, "Call Price":call_prices, "Put Price":put_prices})
print(df.to_markdown(index=False))
