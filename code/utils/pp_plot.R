
# Ploting prior vs. posterior. 
# The inputs are assumed to be the resulting samples from Nimble's MCMC. 
# I expect the user to give equal names to
# hiperparameters on prior and posterior NimbleCode.
# Make sure you have the same monitors as well.
# Remember to load ggplot2.


pp.plot <- function(prior.samples, posterior.samples){
  
  prior <- prior.samples[[1]]
  posterior <- posterior.samples[[1]]
  
  N <- length(prior.samples)
  nomes <- colnames(prior)
  
  if(N >1){
    for(i in 2:N){
      prior <- rbind(prior, prior.samples[[i]])
      posterior <- rbind(posterior, posterior.samples[[i]])
    }
  }
  
  
  #df <- data.frame(prior = prior, posterior = posterior)
  len <- length(prior[,1])

  for(i in 1:length(nomes)){
    df1 <- data.frame(value = prior[,i], type = rep("prior",len))
    df2 <- data.frame(value = posterior[,i], type = rep("posterior",len))
    df <- rbind(df1,df2)
    
    p <- ggplot(data = df, aes(x = value, fill = type,
                               colour = type))+
      geom_density(alpha = 0.4)+
      labs(title = nomes[i])+
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0))+
      theme_bw(base_size = 20)
    
    print(p)
  }
  
}

pp.plot2 <- function(prior, posterior){
  
  len <- length(prior)
  df1 <- data.frame(value = prior, type = rep("prior",len))
  df2 <- data.frame(value = posterior, type = rep("posterior",len))
  df <- rbind(df1,df2)
    
  p <- ggplot(data = df, aes(x = value, fill = type,
                               colour = type))+
    geom_density(alpha = 0.4)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_bw(base_size = 20)
    
  print(p)
}

teste <- function(){
  df <- data.frame(x= 1:5, y= 5:9)
  
  for(i in 1:3){
    p <- ggplot(data = df, aes(x=x, y=y))+geom_line()
    print(p)
  }
}
teste()
