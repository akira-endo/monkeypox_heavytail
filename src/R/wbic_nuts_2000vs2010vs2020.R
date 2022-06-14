## modify data into a form suitable for Stan
x_contact = seq(1,500,1)
N=length(x_contact)
parameters =  c("alpha","beta") # c("llk") 

y_case00 = (rowdf_freq[1,3:(N+2)]) %>% as.numeric() # (rowdf_freq[1,3:(N+2)] + rowdf_freq[5,3:(N+2)])
y_case10 = rowdf_freq[5,3:(N+2)] %>% as.numeric()
y_case20 = rowdf_freq[9,3:(N+2)] %>% as.numeric()

data = list(N=N, x_contact=x_contact, y_case00=y_case00, y_case10=y_case10, y_case20=y_case20)
# specify parameters to monitor

#nuts_fit1 = stan(model_code=Model1,data=data,pars=parameters,iter=3000,thin=1,warmup=2000,chain=5)



## modify data into a form suitable for Stan
y_case00 = (rowdf_freq[3,3:(N+2)]) %>% as.numeric() # (rowdf_freq[1,3:(N+2)] + rowdf_freq[5,3:(N+2)])
y_case10 = rowdf_freq[7,3:(N+2)] %>% as.numeric()
y_case20 = rowdf_freq[11,3:(N+2)] %>% as.numeric()
data = list(N=N, x_contact=x_contact, y_case00=y_case00, y_case10=y_case10, y_case20=y_case20)
# specify parameters to monitor

#nuts_fit2 = stan(model_code=Model2,data=data,pars=parameters,iter=3000,thin=1,warmup=2000,chain=5)



## modify data into a form suitable for Stan
y_case00 = (rowdf_freq[2,3:(N+2)]) %>% as.numeric() # (rowdf_freq[1,3:(N+2)] + rowdf_freq[5,3:(N+2)])
y_case10 = rowdf_freq[6,3:(N+2)] %>% as.numeric()
y_case20 = rowdf_freq[10,3:(N+2)] %>% as.numeric()
data = list(N=N, x_contact=x_contact, y_case00=y_case00, y_case10=y_case10, y_case20=y_case20)
# specify parameters to monitor
#nuts_fit3 = stan(model_code=Model2,data=data,pars=parameters,iter=3000,thin=1,warmup=2000,chain=5)



## modify data into a form suitable for Stan
y_case00 = (rowdf_freq[4,3:(N+2)]) %>% as.numeric() # (rowdf_freq[1,3:(N+2)] + rowdf_freq[5,3:(N+2)])
y_case10 = rowdf_freq[8,3:(N+2)] %>% as.numeric()
y_case20 = rowdf_freq[12,3:(N+2)] %>% as.numeric()
data = list(N=N, x_contact=x_contact, y_case00=y_case00, y_case10=y_case10, y_case20=y_case20)
# specify parameters to monitor
nuts_fit4 = stan(model_code=Model2,data=data,pars=parameters,iter=3000,thin=1,warmup=2000,chain=5)