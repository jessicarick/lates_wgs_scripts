## calculating generation time
## for lates species
## from fishbase estimates
## nov 2021, j. rick

## estimates for l. niloticus from fishbase
Linf = 180 #cm
K = 0.2
M = 0.42

#log (-t0) = -0.3922 - 0.2752 log Linf - 1.038 log K
lt0 = -0.3922 - 0.2752*log(Linf) - 1.038*log(K)
t0 = -exp(lt0)

#Lopt = Linf * (3 / (3 + M/K)) (Beverton 1992)
Lopt = Linf * (3 / (3 + M/K))

#tg = topt = t0 - ln(1 - Lopt / Linf) / K
topt = t0 - log(1-Lopt / Linf) / K
tg1 = topt

tg1 # generation time = 3.58

#------------------------------------------#
## estimates for l. calcarifer from fishbase
Linf = 143.0 #cm
K = 0.1
T = 21.5
M = 10^(0.566 - 0.718 * log(Linf)) + 0.02 * T
        
#log (-t0) = -0.3922 - 0.2752 log Linf - 1.038 log K
lt0 = -0.3922 - 0.2752*log(Linf) - 1.038*log(K)
t0 = -exp(lt0)

#Lopt = Linf * (3 / (3 + M/K)) (Beverton 1992)
Lopt = Linf * (3 / (3 + M/K))

#tg = topt = t0 - ln(1 - Lopt / Linf) / K
topt = t0 - log(1-Lopt / Linf) / K
tg2 = topt

tg2 # generation time = 3.40

#------------------------------------------#
## estimates for l. angustifrons from fishbase
Linf = 195.0 # cm TL
K = 0.25
T = 25
M = 10^(0.566 - 0.718 * log(Linf)) + 0.02 * T

#log (-t0) = -0.3922 - 0.2752 log Linf - 1.038 log K
lt0 = -0.3922 - 0.2752*log(Linf) - 1.038*log(K)
t0 = -exp(lt0)

#Lopt = Linf * (3 / (3 + M/K)) (Beverton 1992)
Lopt = Linf * (3 / (3 + M/K))

#tg = topt = t0 - ln(1 - Lopt / Linf) / K
topt = t0 - log(1-Lopt / Linf) / K
tg3 = topt

tg3 # generation time = 2.99

#------------------------------------------#
## estimates for l. microlepis from fishbase
Linf = 108.0 # cm TL
K = 0.34
T = 25
M = 10^(0.566 - 0.718 * log(Linf)) + 0.02 * T

#log (-t0) = -0.3922 - 0.2752 log Linf - 1.038 log K
lt0 = -0.3922 - 0.2752*log(Linf) - 1.038*log(K)
t0 = -exp(lt0)

#Lopt = Linf * (3 / (3 + M/K)) (Beverton 1992)
Lopt = Linf * (3 / (3 + M/K))

#tg = topt = t0 - ln(1 - Lopt / Linf) / K
topt = t0 - log(1-Lopt / Linf) / K
tg4 = topt

tg4 # generation time = 2.69

#------------------------------------------#
## estimates for l. mariae from fishbase
Linf = 78.0 # cm TL
K = mean(0.253,0.120,0.260,0.250)
T = 26
M = 10^(0.566 - 0.718 * log(Linf)) + 0.02 * T

#log (-t0) = -0.3922 - 0.2752 log Linf - 1.038 log K
lt0 = -0.3922 - 0.2752*log(Linf) - 1.038*log(K)
t0 = -exp(lt0)

#Lopt = Linf * (3 / (3 + M/K)) (Beverton 1992)
Lopt = Linf * (3 / (3 + M/K))

#tg = topt = t0 - ln(1 - Lopt / Linf) / K
topt = t0 - log(1-Lopt / Linf) / K
tg5 = topt

tg5 # generation time = 2.70

#------------------------------------------#
## estimates for l. stappersii from fishbase
Linf = 50.5 # cm TL
K = mean(0.398,0.4,0.4,0.39,0.39,0.38,0.25,0.43,0.35,0.4,0.4,0.355,0.36)
T = 26
M = 10^(0.566 - 0.718 * log(Linf)) + 0.02 * T

#log (-t0) = -0.3922 - 0.2752 log Linf - 1.038 log K
lt0 = -0.3922 - 0.2752*log(Linf) - 1.038*log(K)
t0 = -exp(lt0)

#Lopt = Linf * (3 / (3 + M/K)) (Beverton 1992)
Lopt = Linf * (3 / (3 + M/K))

#tg = topt = t0 - ln(1 - Lopt / Linf) / K
topt = t0 - log(1-Lopt / Linf) / K
tg6 = topt

tg6 # generation time = 2.38


#------------------------------------------#
mean(c(tg1,tg2,tg3,tg4,tg5,tg6)) # overall average 2.96
range(c(tg1,tg2,tg3,tg4,tg5,tg6))
