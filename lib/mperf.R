mperf     = function(sim,obs,vnam,dchart,outidx){
  
  #--------------------------------------------------#
  #------------- Performance function ---------------#
  #--- Compute statistical indexes of performance ---#
  #--------------------------------------------------#
  #
  # Decription:
  # sim     - Simulated values          [Real]
  # obs     - Observed values           [Real]
  # vnam    - Name of variable          [String]
  # dchart  - Display Chart?            [T or F]
  # outidx  - Output peformance indexes [List]
  # 
  # Murilo Vianna
  # source: https://github.com/Murilodsv/R-scripts/blob/master/mperf.r
  #
  # Literature: 
  # Brun, F., Wallach, D., Makowski, D., & Jones, J. W. (2006). 
  # Working with dynamic crop models: evaluation, analysis, 
  # parameterization, and applications. Elsevier.
  #--------------------------------------------------#
  
  if(missing(sim)){stop("Missing sim argument")}
  if(missing(obs)){stop("Missing obs argument")}
  if(missing(dchart)){dchart = T}
  if(missing(outidx)){outidx = c("bias","mse","rmse","mae","rrmse","rmae","ef","r","r2","d","cc","a","b","mi.sim","sd.sim","cv.sim","mi.obs","sd.obs","cv.obs","n")}
  if(missing(vnam)){
    warning("Missing vnam argument: vnam set to none.")
    vnam = ""
  }
  
  #--- all outputs
  if(outidx[1] == "all"){outidx = c("bias","mse","rmse","mae","rrmse","rmae","ef","r","r2","d","cc","a","b","mi.sim","sd.sim","cv.sim","mi.obs","sd.obs","cv.obs","n")}
  
  #--- Check Input data
  sim = as.numeric(sim)
  obs = as.numeric(obs)
  
  if(length(sim) != length(obs)){stop("Vector length of Simulated and Observed do not match.")}
  if(length(sim[is.na(sim)]) > 0 |  length(obs[is.na(obs)]) > 0){
    
    #--- remove NA values
    df_so = na.omit(data.frame(sim,obs))
    sim = df_so$sim
    obs = df_so$obs
    warning("NA values ignored in Simulated or Observed data")
  }
  
  #--- Statistical indexes
  n     = length(obs)
  mi.sim= mean(sim)
  mi.obs= mean(obs)
  fit   = lm(sim~obs)
  bias  = (1/n) * sum(sim-obs)
  mse   = (1/n) * sum((sim-obs)^2)
  rmse  = sqrt(mse)
  mae   = (1/n) * sum(abs(sim-obs))
  rrmse = rmse / mean(obs)
  rmae  = (1/length(obs[obs>0])) * sum(abs(sim[obs>0]-obs[obs>0])/abs(obs[obs>0]))
  ef    = 1 - (sum((sim-obs)^2) / sum((obs-mean(obs))^2))
  r     = sum((obs-mean(obs))*(sim-mean(sim)))/sqrt(sum((obs-mean(obs))^2)*sum((sim-mean(sim))^2))
  r2    = r^2
  d     = 1 - (sum((sim-obs)^2) / sum((abs(sim-mean(obs))+abs(obs-mean(obs)))^2))
  if(length(unique(sim)) > 1){
    a     = fit$coefficients[1]
    b     = fit$coefficients[2]
    if(any(is.na(c(a,b)))){
      outidx = outidx[outidx != "a" & outidx != "b"]
    }
  }else{
    outidx = outidx[outidx != "a" & outidx != "b"]
    a = NA
    b = NA
  }
  
  sigma.obs.sim = (1 / n) * sum((obs - mi.obs) * (sim - mi.sim))
  sigma.obs.2   = (1 / n) * sum((obs - mi.obs) ^ 2)
  sigma.sim.2   = (1 / n) * sum((sim - mi.sim) ^ 2)
  cc    = (2 * sigma.obs.sim) / (sigma.obs.2 + sigma.sim.2 + (mi.sim - mi.obs)^2)
  
  sd.obs = sd(obs)
  sd.sim = sd(sim)
  cv.obs = sd.obs / mi.obs
  cv.sim = sd.sim / mi.sim
  
  if(dchart){
    #--- Chart Sim ~ Obs
    varlab = vnam 
    
    mindt = min(obs,sim)
    maxdt = max(obs,sim)
    #--- Ploting limits 
    pllim = c(mindt-0.1*(maxdt-mindt),maxdt+0.1*(maxdt-mindt))
    xx = seq(min(obs),max(obs),length = (max(obs)-min(obs))*1000)
    z = summary(fit)
    
    plot(sim~obs,
         ylab = paste("Sim - ",varlab,sep = ""),
         xlab = paste("Obs - ",varlab,sep = ""),
         ylim = pllim,
         xlim = pllim)
    
    lines(xx, predict(fit, data.frame(obs=xx)),
          col = "black",
          lty = 1,
          lwd = 1.5)
    
    l11 = seq(pllim[1]-0.5*(maxdt-mindt), pllim[2] + 0.5 * (maxdt-mindt),length = 1000)
    
    lines(l11*1~l11,
          col = "red",
          lty = 2,
          lwd = 1.5)
  }
  
  perf = data.frame(model = vnam,
                    bias,
                    mse,
                    rmse,
                    mae,
                    rrmse,
                    rmae,
                    ef,
                    r,
                    r2,
                    d,
                    cc,
                    a,
                    b,
                    mi.sim,
                    sd.sim,
                    cv.sim,
                    mi.obs,
                    sd.obs,
                    cv.obs,
                    n)
  
  rownames(perf) = c()
  
  return(perf[,c("model",outidx)])
  
}
