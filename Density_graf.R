# Graphics: Probability Density Function EGEDD ------------------------------- #

library(tidyverse)
library(latex2exp)
library(ggplot2)
library(tibble)
library(purrr)

# ---------------------------------------------------------------------------- #

#Probability Density Function UMW1----
d_UMW1<-function(y,alpha,gamma,lambda)
{
  f<- (alpha/(log(y)*y^(lambda+1)))*((-log(y))^gamma)*(lambda*log(y)-gamma)*exp(-alpha*((-log(y))^gamma)*y^(-lambda))
  return(f)
}
#Cumulative Function UMW1----
p_UMW1<-function(y,alpha,gamma,lambda)
{
  f<- exp(-alpha*((-log(y))^gamma)*(y^(-lambda)))
  return(f)
}

#Probability Density Function UMW1----
d_UMW_q<-function(y,mu,gamma,lambda,tau=0.5)
{
  f<- (-mu^lambda * log(tau) * (-log(y))^gamma * (lambda * log(y) - gamma)) / (y^(lambda + 1) * (-log(mu))^gamma * log(y)) * exp(mu^lambda * log(tau) * (-log(y))^gamma / (y^lambda * (-log(mu))^gamma))
  return(f)
}



# ---------------------------------------------------------------------------- #

graphic_density<-function(par,x,data,legend_g=r'($\alpha = %f$)',lim_x,y_max,legend_pos=c(0.9, 0.8),base_size=12){
  #http://sape.inf.usi.ch/quick-reference/ggplot2/linetype
  #myDot<-c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "1F", "F1", "4C88C488", "12345678")
  myDot<-c("solid", "longdash","twodash","dashed", "dotdash", "dotted")
  my_color<-rep("black",6)
  #my_color<-viridis::cividis(15)
  last_plot<-ggplot2::ggplot(data,aes(x = x, y = v, colour =as.factor({{par}}), linetype=as.factor({{par}})))+
    xlim(lim_x)+ylim(0,y_max)+geom_line(size = 0.8)+
    scale_linetype_manual(name= "", values = myDot,labels=lapply(sprintf(legend_g, {{par}}), TeX))+
    scale_color_manual(name="",values=my_color,labels=lapply(sprintf(legend_g, {{par}}), TeX))+
    labs(x="y",y = "Probability Density Function",title="",linetype = "")+theme_bw(base_size = base_size)+
    guides(linetype = guide_legend(ncol = 1)) +
    theme(legend.position = legend_pos,panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),legend.background = element_blank(),
          legend.key = element_blank())
  last_plot
}

# UMW1 =========================================================================

alpha<-c(0.1,0.2,0.5,1.5,3.5,5.5)
gamma<-c(0.1,0.5,0.8,1.5,3.5,5.5)
lambda<-c(0.1,0.3,0.7,1.5,3.5,5.5)
x<-seq(.0001,.9999,length.out = 1000)

# Alpha ---------------------------------------------------------------------- #

dataa<-purrr::map_df(alpha, ~ tibble::tibble(v=d_UMW1(x,.,gamma=1.2,lambda=0.7), x=x, alpha=.))
last_plota<-graphic_density(par = alpha,x = x,data = dataa,legend_g = r'($\alpha = %f$)',
                            lim_x = c(0,1),y_max = 4.8,legend_pos = c(0.45, 0.77),base_size = 13)
last_plota
ggsave("pdf_UMW1_valpha.pdf",plot = last_plota, width = 11.5, height = 9, units = "cm")

# gamma ---------------------------------------------------------------------- #

datab<-purrr::map_df(gamma, ~ tibble::tibble(v=d_UMW1(x,.,alpha=0.4, lambda=0.3), x=x, gamma=.))
last_plotb<-graphic_density(par = gamma,x = x,data = datab,legend_g = r'($\gamma = %f$)',
                            lim_x = c(0,1),y_max = 6.5,legend_pos = c(0.7, 0.77),base_size = 13)
last_plotb
ggsave("pdf_UMW1_vgamma.pdf",plot = last_plotb, width = 11.5, height = 9, units = "cm")

# lambda --------------------------------------------------------------------- #

datac<-purrr::map_df(lambda, ~ tibble::tibble(v=d_UMW1(x,.,alpha=0.6,gamma=1.5), x=x, lambda=.))
last_plotc<-graphic_density(par = lambda,x = x,data = datac,legend_g = r'($\lambda = %f$)',
                            lim_x = c(0,1),y_max = 5.5,legend_pos = c(0.15, 0.77),base_size = 13)
last_plotc
ggsave("pdf_UMW1_vlambda.pdf",plot = last_plotc, width = 11.5, height = 9, units = "cm")


# UMW1 q =======================================================================

mu<-c(0.1,0.3,0.5,0.7,0.9)
gamma<-c(0.1,0.5,1.5,3.5,5.5)
lambda<-c(0.1,0.5,1.5,3.5,5.5)
x<-seq(.0001,.9999,length.out = 1000)

# Alpha ---------------------------------------------------------------------- #

dataa<-purrr::map_df(mu, ~ tibble::tibble(v=d_UMW_q(x,.,gamma=1.5,lambda=0.4), x=x, mu=.))
last_plota<-graphic_density(par = mu,x = x,data = dataa,legend_g = r'($\mu = %f$)',
                            lim_x = c(0,1),y_max = 7.2,legend_pos = c(0.25, 0.8))
last_plota
ggsave("pdf_UMW1q_vmu.pdf",plot = last_plota, width = 12, height = 9, units = "cm")

# gamma ---------------------------------------------------------------------- #

datab<-purrr::map_df(gamma, ~ tibble::tibble(v=d_UMW_q(x,.,mu=0.5, lambda=0.3), x=x, gamma=.))
last_plotb<-graphic_density(par = gamma,x = x,data = datab,legend_g = r'($\gamma = %f$)',
                            lim_x = c(0,1),y_max = 6.5,legend_pos = c(0.8, 0.8))
last_plotb
ggsave("pdf_UMW1q_vgamma.pdf",plot = last_plotb, width = 12, height = 9, units = "cm")

# lambda --------------------------------------------------------------------- #

datac<-purrr::map_df(lambda, ~ tibble::tibble(v=d_UMW_q(x,.,mu=0.4,gamma=1.2), x=x, lambda=.))
last_plotc<-graphic_density(par = lambda,x = x,data = datac,legend_g = r'($\lambda = %f$)',
                            lim_x = c(0,1),y_max = 6.7,legend_pos = c(0.88, 0.8))
last_plotc
ggsave("pdf_UMW1q_vlambda.pdf",plot = last_plotc, width = 12, height = 9, units = "cm")

