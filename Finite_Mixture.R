
packages = c("nimble","codetools","devtools","Rcpp","dplyr",
             "ggplot2","parallel","doParallel","doSNOW","ggmcmc","reshape2","job")
invisible(xfun::pkg_attach(packages))


#defining required distribution
## 1 Generalized Pareto Distribution


#Density Function
dgpd <- nimbleFunction(
  run = function(x = double(0),u = double(0),
                 sigma = double(0),xi = double(0),
                 log = integer(0, default = 0))
  {
    returnType(double(0))
    z = (x - u) / sigma
    if (x < u)
      logProb = -Inf
    else{
      if (xi == 0)
        logProb = -z - log(sigma)
      else
      {
        if (xi < 0)
        {
          if (x > u - (sigma / xi))
            logProb = -Inf
          else
            logProb = -(1 + (1 / xi)) * log(1 + (xi * z)) - log(sigma)
        }
        else
          logProb = -(1 + (1 / xi)) * log(1 + (xi * z)) - log(sigma)
      }
    }
    if (log)
      return(logProb)
    else
      return(exp(logProb))
  }
)
compileNimble(dgpd)
assign('dgpd', dgpd, .GlobalEnv)

#Generating Function
rgpd <- nimbleFunction(
  run=function(n = integer(0),u=double(0)
               ,sigma=double(0),xi=double(0))
  {
    if(n != 1) print("rgpdonly allows n = 1;")
    else{
      y=runif(1)
      if(xi==0) x=u-sigma*log(y)
      else x=u+((sigma/xi)*((y^(-xi))-1))
    }
    return(x)
    returnType(double(0))
  })
compileNimble(rgpd)
assign('rgpd', rgpd, .GlobalEnv)




#Cumulative sum
cum_sum=nimbleFunction(run=function(p=double(1))
{
  returnType(double(1))
  for(j in 2:length(p))
    p[j]=p[j-1]+p[j]
  return(p)
})
compileNimble(cum_sum)
assign('cum_sum',cum_sum,.GlobalEnv)


#Density Function
dTotal<- nimbleFunction(
  run = function(x=double(0),alpha=double(0),lambda=double(1)
                 ,mu=double(1),u=double(0),sigma=double(0)
                 ,xi=double(0),log = integer(0, default = 0)) 
  {
    returnType(double(0))
    d=length(lambda)
    w=rdirch(1,rep(alpha/d,d))
    if(x<u){
      p1=dgamma(x,lambda[1:d],mu[1:d])*w[1:d]
      prob=sum(p1)
    }
    else
    {
      p2=pgamma(u,lambda[1:d],mu[1:d])*w[1:d]
      P2=sum(p2)
      prob=(1-P2)*dgpd(x,u,sigma,xi)
    }
    if(!log) return(prob)
    else return(log(prob))
  })
compileNimble(dTotal) 
assign('dTotal', dTotal, .GlobalEnv)

#Generating Function 
rTotal <- nimbleFunction(run=function(n = integer(0)
                                      ,alpha=double(0),lambda=double(1)
                                      ,mu=double(1),u=double(0)
                                      ,sigma=double(0),xi=double(0))
{
  returnType(double(0))
  if(n != 1) print("rTotal only allows n = 1;")
  else{
    y=runif(1)
    d=length(lambda)
    w=rdirch(1,rep(alpha/d,d))
    W=cum_sum(w)
    g=pgamma(u,lambda,mu)
    p=g*w
    P=sum(p)
    for(i in 2:d)
    {
      k=runif(1)
      if(W[i-1]<k & W[i]>k) x=rgamma(1,lambda[i],mu[i])
      else x=rgamma(1,lambda[1],mu[1])
    }
    a=(1-y)/(1-P)
    if(xi==0) r=u-sigma*log(a)
    else r=u+((sigma/xi)*(a^(-xi)-1))
    if(y<P) return(x)
    else return(r)
  }
  
})   

compileNimble(rTotal)
assign('rTotal', rTotal, .GlobalEnv)

# Function to Generate Data(r=#replication,N=#data points,p=tail_index)
data_sim <- function(r=10,N=200,p=0.9,lambda=c(10,6),mu=c(4,0.7)
                     ,u=11,sigma=3,xi=0.4)
{
  Y = data.frame(matrix(ncol=r,nrow=N))
  n_rep=1
  repeat{
    for(i in 1:N)
    {
      if(runif(1)<0.9)
        Y[i,n_rep]=ifelse(runif(1)<0.5,rgamma(1,lambda[1],mu[1])
                          ,rgamma(1,lambda[2],mu[2]))
      else
        Y[i,n_rep]=ifelse(runif(1)<0.5,rgamma(1,lambda[1],mu[1])
                          +rgpd(1,u,sigma,xi)
                          ,rgamma(1,lambda[2],mu[2]))
      +rgpd(1,u,sigma,xi)
    }
    n_rep=n_rep+1
    if(n_rep>r){break}
  }
  names(Y) <- NULL  
  return(Y)
}
assign('data_sim',data_sim,.GlobalEnv)

# function to generate true density
dTrue <- function(x,lambda=c(10,6),mu=c(4,0.7)
                  ,u=11,sigma=3,xi=0.4)
{
  if(x<u)
    density <- (dgamma(x,lambda[1],mu[1])+dgamma(x,lambda[2],mu[2]))/2
  else 
    density = (1-(pgamma(u,lambda[1],mu[1])+pgamma(u,lambda[2],mu[2]))/2)*dgpd(x,u,sigma,xi)
  return(density)
}
assign("dTrue",dTrue,.GlobalEnv)

# Function to generate sample from True density
rTrue <- function(n,lambda=c(10,6),mu=c(4,0.7)
                  ,u=11,sigma=3,xi=0.4)
{
  r=array()
  for(i in 1:n)
  {
    if(runif(1)<0.9)
      r[i]=ifelse(runif(1)<0.5,rgamma(1,lambda[1],mu[1])
                        ,rgamma(1,lambda[2],mu[2]))
    else
      r[i]=ifelse(runif(1)<0.5,rgamma(1,lambda[1],mu[1])
                        +rgpd(1,u,sigma,xi)
                        ,rgamma(1,lambda[2],mu[2]))+rgpd(1,u,sigma,xi)
  }
  return(r)
}
assign("rTrue",rTrue,.GlobalEnv)

# Generating data
Y=data_sim()
#Estimating Prior sd for U
sigma_prior = function(y,sigma)
{
  root = function(sigma){
    a1=quantile(y,0.99)
    a2=quantile(y,0.5)
    fun <- pnorm(a1,11,sigma) -pnorm(a2,11,sigma)-0.99
    return(fun)
  }
  uniroot(root,c(0,100))$root
}
sigma_u <- mean(sapply(Y,sigma_prior));sigma_u
max_y = max(Y)

#model Code
code=nimbleCode(
{
    #parameters for Central part and Tail part (DPMG)
    alpha ~ dgamma(1,1)
    for(i in 1:m)
    {
      lambda[i] ~ dgamma(0.001,0.001)
      mu[i] ~ dgamma(0.001,0.001)
    }
    #additional Parameters for Tail Part
    u ~ dnorm(10,sd=1)
    xi ~ T(dnorm(0.4,sd=0.5),-0.5,Inf)
    sigma ~ dgamma(1,1)
    #Distribution for data
    for(i in 1:n)
    {
      y[i] ~ dTotal(alpha,lambda[1:m],mu[1:m],u,sigma,xi)
    }
  })

#MCMC simulation

Sim <- function(y,n_mix=5,n_iter=5000,n_burnin=2500,n_chains=1,n_thin=1)
{
  data = list(y=y)
  ##Setting Up model
  # Model Constants
  consts=list(n=length(data$y),m=n_mix)
  # Parameter initialization,
  inits=list(lambda=rexp(consts$m, 1),mu=rexp(consts$m, 1)
             ,alpha = 1,xi=0.4,sigma=3,u=quantile(data$y,0.9))
  t1=Sys.time()
  nimble_out <- nimbleMCMC(code, constants = consts, data = data,
                           monitors = c("alpha","lambda" ,"mu","u","sigma","xi"),
                           inits = inits,niter = n_iter,nburnin = n_burnin,
                           nchains=n_chains,thin = n_thin,progressBar = FALSE,
                           samplesAsCodaMCMC = TRUE,summary = TRUE)
  t2=Sys.time()
  if(n_chains==1)
  { 
    estimates = nimble_out$summary[,1]
    summary = nimble_out$summary[c("alpha","sigma","u","xi"),]
    summary_all = nimble_out$summary[,1]
  }
  else
  {
    estimates = nimble_out$summary$all.chains[,1]
    summary = nimble_out$summary$all.chains[c("alpha","sigma","u","xi"),]
    summary_all = nimble_out$summary$all.chains[,1]
  }
  
  t=t2-t1
  #Storing Result
  my.list = list()
  my.list$run_time = as.numeric(t)
  my.list$model_summary = summary
  # estimated parameters
  alpha = estimates["alpha"]
  
  lambda = summary_all[seq(2,n_mix+1)]
  my.list$L=lambda
  mu = summary_all[seq(n_mix+2,(2*n_mix)+1)]
  my.list$M = mu
  u = estimates["u"]
  sigma = estimates["sigma"]
  xi = estimates["xi"]
  # MCMC plot
  S=ggs(nimble_out$samples[,c('alpha','u','sigma','xi')])
  my.list$trace_plot = ggs_traceplot(S,greek = TRUE)
  my.list$density_plot = ggs_density(S,greek = TRUE)
  #calculating Integrated Standard error
  grid = seq(min(y),max(y), by = 0.01)
  f_hat = sapply(grid, dTotal,alpha,lambda,mu,u,sigma,xi)
  f_true = sapply(grid, dTrue)
  ISE = mean((f_hat-f_true)^2)
  my.list$Integrated_Squared_error = ISE
  # density comparison
  grid1 = seq(min(y),max(y), length.out = length(y))
  f_hat = sapply(grid1, dTotal,alpha,lambda,mu,u,sigma,xi)
  f_true = sapply(grid1, dTrue)
  output = data.frame(grid1,f_true,f_hat)
  names(output) = c("grid","True","Estimated")
  p = ggplot(output,aes(x=grid))+
    geom_line(aes(y=True,color="True"))+
    geom_line(aes(y=Estimated,color="Posterior"))+
    xlab('Grid')+ylab('Density') 
  my.list$plot = p
  my.list$y = y
  return(my.list)
}
assign("Sim",Sim,.GlobalEnv)
#v=sapply(Y[,c(2,3,4)],Sim,n_chains=1,n_mix = 5,n_thin = 1,n_iter=300,n_burnin=150)
#out=v
# unlist(v['y',])
# Result

Analyze=function(out,n_mix=5,u=11,sigma=3,xi=0.4,
                 n_iter=1000,n_burnin=0,n_chains=1,n_thin=1)
{
  # estimated parameter value
  a_m=array()
  U_m=array()
  x_m=array()
  s_m=array()
  ise=array()
  diagnostics=list()
  n_rep=ncol(out)
  for(i in 1:n_rep)
  {
    a_m[i] = out['model_summary',i]$model_summary["alpha","Mean"]
    U_m[i] = out['model_summary',i]$model_summary["u","Mean"]
    s_m[i] = out['model_summary',i]$model_summary["sigma","Mean"]
    x_m[i] = out['model_summary',i]$model_summary["xi","Mean"]
    ise[i] = out['Integrated_Squared_error',i]
    diagnostics$model_summary[[i]]= out['model_summary',i]
    diagnostics$trace_plots[[i]] = out['trace_plot',i]
    diagnostics$density_plots[[i]] = out['density_plot',i]
    if(n_chains>1)
    {
      diagnostics$R_hat.plots[[i]] = out['R_hat',i]
      diagnostics$geweke.plots[[i]] = out['geweke',i]
    }
    diagnostics$true_vs_posterior[[i]] = out['plot',i]
    diagnostics$run_time[[i]] = out['run_time',i]
  }
  # parameter Estimates
  alpha_hat=mean(a_m)
  lambda_hat=colMeans(matrix(as.numeric(unlist(out['L',])),ncol = n_mix,byrow=T))
  mu_hat=colMeans(matrix(as.numeric(unlist(out['M',])),ncol = n_mix,byrow=T))
  u_hat=mean(U_m)
  sigma_hat=mean(s_m)
  xi_hat=mean(x_m)
  
  # Calculating Bias
  bias_u = u_hat-u
  bias_sigma = sigma_hat-sigma
  bias_xi = xi_hat-xi
  # calculating empirical SE
  empirical_SE_u = sd(a_m)
  empirical_SE_sigma = sd(s_m)
  empirical_SE_xi = sd(x_m)
  
  # number of effective sample
  n_total = (n_iter-n_burnin)*(n_chains/n_thin)
  # calculating Monte-Carlo standard error
  MC_SE_u = empirical_SE_u/sqrt(n_total)
  MC_SE_sigma = empirical_SE_sigma/sqrt(n_total)
  MC_SE_xi = empirical_SE_xi/sqrt(n_total)
  
  #storing result
  final = list()
  final$MISE = mean(as.numeric(ise))
  final$total.samples.used = n_total
  final$summary.mcmc = matrix(c(u_hat,sigma_hat,xi_hat,bias_u,bias_sigma,bias_xi
                                ,empirical_SE_u,empirical_SE_sigma
                                ,empirical_SE_xi,MC_SE_u,MC_SE_sigma,MC_SE_xi),nrow=3)
  colnames(final$summary.mcmc) = c("Estimates","Bias","Empirical_SE","MC_SE")
  rownames(final$summary.mcmc) = c("u","sigma","xi")
  final$mcmc.diagnostics = diagnostics
  final$posterior_estimate = list(alpha_hat,lambda_hat,mu_hat,u_hat,sigma_hat,xi_hat)
  names(final$posterior_estimate ) = c("alpha","lambda","mu","u","sigma","xi")
  # density Plot of parameters
  samples = data.frame('alpha'=a_m,'u'=U_m,'sigma'= s_m,'Xi'=x_m)
  posterior = list()
  posterior$alpha = ggplot(samples,aes(x=alpha),after_stat(density))+
                    geom_histogram(aes(y=after_stat(density)),color="blue",fill="lightgreen")+
                      geom_density(color="red",linewidth=1)
                      
  posterior$u = ggplot(samples,aes(x=u))+
                   geom_histogram(aes(y=after_stat(density)),color="blue",fill="lightgreen")+
                      geom_density(color="red",linewidth=1)
  posterior$sigma = ggplot(samples,aes(x=sigma),after_stat(density))+
                      geom_histogram(aes(y=after_stat(density)),color="blue",fill="lightgreen")+
                      geom_density(color="red",linewidth=1)
                      
  posterior$xi=  ggplot(samples,aes(x=Xi),after_stat(density))+
                 geom_histogram(aes(y=after_stat(density)),color="blue",fill="lightgreen")+
                      geom_density(color="red",linewidth=1)
                      
  final$Posterior_parameter_density = posterior

  return(final)
}

assign('Analyze',Analyze,.GlobalEnv)
 # output = Analyze(v)
 # output$Posterior_parameter_density$u
# rm(Report)

Report = function(n_rep=10,N=100,lambda=c(10,6),mu=c(4,0.7),p=0.9,n_mix=4,
                  u=11,sigma=3,xi=0.4,n_iter=100,n_burnin=0,n_chains=1,
                  n_thin=1,seed=10,n_cores=11)
{
  
  set.seed(seed)
  t1=Sys.time()
  Y=data_sim(n_rep,N,p,lambda,mu,u,sigma,xi)
  
  ## parallel computing Using Apply
  my.cluster <- makeCluster(n_cores)
  registerDoParallel(my.cluster)
  
  clusterEvalQ(my.cluster, 
               {library(nimble)
                 library(codetools)
                 library(devtools)
                 library(Rcpp)
                 library(dplyr)
                 library(ggplot2)
                 library(parallel)
                 library(doParallel)
                 library(foreach)
                 library(doSNOW)
                 library(ggmcmc)
                 library(reshape2)
               })
  clusterExport(my.cluster,c("dgpd","rgpd","dTrue","rTrue","cum_sum","dTotal","rTotal","code"
                             ),envir = .GlobalEnv)
  out <- parSapply(my.cluster,Y,Sim,n_mix,n_iter,n_burnin,n_chains,n_thin)
  stopCluster(my.cluster)
  output = Analyze(out,n_mix,u,sigma,xi,n_iter,n_burnin,n_chains,n_thin)
  t=Sys.time()-t1
  my.list = list()
  my.list$Mean_Int_Sq_error = output$MISE
  my.list$Total_execution_time = t
  my.list$model_Summary = output$summary.mcmc
  my.list$mcmc.diagnostics = output$mcmc.diagnostics
  my.list$Posterior_density = output$Posterior_parameter_density
  # density comparison
  my.list$lambda_hat = output$posterior_estimate$lambda
  my.list$mu_hat = output$posterior_estimate$mu
  size = output$total.samples.used
  grid = seq(0,40, length.out = size)
  f_hat = as.numeric(sapply(grid, dTotal,output$posterior_estimate$alpha,
                 output$posterior_estimate$lambda,output$posterior_estimate$mu,
                 output$posterior_estimate$u,output$posterior_estimate$sigma,
                 output$posterior_estimate$xi))
  
  f_true = as.numeric(sapply(grid, dTrue))
  r_true = rTrue(size)

  p = ggplot()+xlim(0,40)+ylim(0,max(f_true)+0.001)+xlab('')+
    geom_histogram(aes(r_true,after_stat(density))
                   ,color="black",fill="beige",bins = 40)+
    geom_line(aes(x=grid,y=f_true ),color="red",linewidth=1)+
    geom_smooth(aes(x=grid,f_hat),color = "blue",linetype="dotted")+
    geom_segment(aes(x = u , y = 0, xend = u, yend = .27),color='red',linewidth=1)+
    geom_segment(aes(x = u_hat , y = 0, xend = u_hat, yend = .27)
                 ,linetype="dotted",color='blue',linewidth=1.2)+ guides(col = "none")+ 
    ylab('Density')+ theme(legend.position = "none")  
    my.list$Density_Plots = p
  
  # plot zoom
  p1 = ggplot()+xlim(8,20)+ylim(0,0.1)+xlab('')+
    geom_histogram(aes(r_true,after_stat(density))
                   ,color="black",fill="beige",bins = 40)+
    geom_line(aes(x=grid,y=f_true ),color="red",linewidth=1)+
    geom_smooth(aes(x=grid,f_hat),color = "blue",linetype="dotted")+
    geom_segment(aes(x = u , y = 0, xend = u, yend = .27),color='red',linewidth=1)+
    geom_segment(aes(x = u_hat , y = 0, xend = u_hat, yend = .27)
                 ,linetype="dotted",color='blue',linewidth=1.2)+ guides(col = "none")
    ylab('Density')+ theme(legend.position = "none")  
  my.list$zoomed_Density_Plots = p1
  
  return(my.list)
}
rm(report1)
report1 = Report(n_rep = 200,N=200,n_iter = 8500,n_burnin=2500,n_chains=1,
                 seed = 10,n_mix=2,n_cores = 9,n_thin = 3)

report1$Mean_Int_Sq_error
report1$Total_execution_time
report1$zoomed_Density_Plots
report1$Density_Plots
report1$model_Summary
report1$Posterior_density
report1$lambda_hat
report1$mu_hat
################################
report2 = Report(n_rep = 200,N=200,n_iter = 8500,n_burnin=2500,n_chains=1,
                 seed = 10,n_mix=2,n_cores = 8,n_thin = 3)

report2 $Mean_Int_Sq_error
report2 $Total_execution_time
report2 $zoomed_Density_Plots
report2 $Density_Plots
report2 $model_Summary
report2 $Posterior_density
######################################################################################
report3  = Report(n_rep = 200,N=200,n_iter = 8500,n_burnin=2500,n_chains=1,
                 seed = 10,n_mix=2,n_cores = 8,n_thin = 3)

report3 $Mean_Int_Sq_error
report3 $Total_execution_time
report3 $zoomed_Density_Plots
report3 $Density_Plots
report3 $model_Summary
report3 $Posterior_density
###################################################################
report4 = Report(n_rep = 200,N=200,n_iter = 8500,n_burnin=2500,n_chains=1,
                 seed = 10,n_mix=15,n_cores = 8,n_thin = 3)

report4 $Mean_Int_Sq_error
report4 $Total_execution_time
report4 $zoomed_Density_Plots
report4 $Density_Plots
report4 $model_Summary
report4 $Posterior_density
#########################################################
report5 = Report(n_rep = 200,N=200,n_iter = 8500,n_burnin=2500,n_chains=1,
                 seed = 10,n_mix=20,n_cores = 8,n_thin = 3)

report5 $Mean_Int_Sq_error
report5 $Total_execution_time
report5 $zoomed_Density_Plots
report5 $Density_Plots
report5 $model_Summary
report5 $Posterior_density

