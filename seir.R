library(deSolve)
library(ggplot2)
library(svglite)
library(devtools)
library(roxygen2)
library(knitr)
library(testthat)
library(rstudioapi)


#' ODE solution
#'
#' This function solves an ODE system based on the given parameter list.
#' @param user created parameters
#' @return solution of an ODE function system
#' @export
model <- function(t, y, param){
  S <- y[1]
  E <- y[2]
  I1 <- y[3]
  I2 <- y[4]
  R <- y[5]
  N <- param["N"]
  beta <- param["beta"]
  mu <- param["mu"]
  gamma <-param["gamma"]
  lamda <-param["lamda"]
  #传染病数学模型
  dSt <- mu * (N - S) - beta * S * I1/N
  dEt <- beta * S * I1/N - mu * E-lamda*E
  dI1t <-   - (mu + gamma) * I1+lamda*E
  dI2t <-   - (mu + gamma) * I1+lamda*E*E
  dRt <- gamma * I1 - mu * R
  #求解结果整合成向量表达
  outcome <- c(dSt, dEt, dI1t, dI2t, dRt)
  #返回常微分方程系统求解结果
  list(outcome)
}
#设置评估参数的初始值
times <- seq(0, 156, by = 1/7)
param <- c(mu = 0.000, lamda = 0.03, beta = 4, gamma=0.1, N = 1)
init <- c(S = 0.9999, E = 0.00008, I1 = 0.00002, I2 = 0.00001, R = 0)
#调用常微分方程求解函数，传入初始条件，评估时间，模型以及参数信息



result <-  deSolve::ode(y=init, times=times, func=model, parms = param)
result <- as.data.frame(result)
tail(round(result, 3.6),10)

#结果画图
#' @export
seirplot <- ggplot2::ggplot(data=result)+
  ggplot2::geom_line(ggplot2::aes(x=time, y=S,col="S"), lwd=2)+
  ggplot2::geom_line(ggplot2::aes(x=time, y=I1,col="I1"), lwd=2)+
  ggplot2::geom_line(ggplot2::aes(x=time, y=I2,col="I2"), lwd=2)+
  ggplot2::geom_line(ggplot2::aes(x=time, y=R,col="R"), lwd=2)+
  ggplot2::geom_line(ggplot2::aes(x=time, y=E,col="E"), lwd=2)+
  ggplot2::labs(x = "Time",y = "Ratio")+
  ggplot2::scale_color_manual(name = "SEIR",
     values = c("S" = "orange", "E" = "purple", "I1" = "red", "I2" = "blue" ,"R" = "green"))
#绘制仿真结果并保存为矢量文件
seirplot
ggplot2::ggsave(seirplot, file="seir.pdf", width=7, height=6)
ggplot2::ggsave(seirplot, file="seir.svg", width=7, height=6)

getwd()
#[1] "D:\中本一郎大三下\数据挖掘课后作业3\seirPackage"
devtools::document()
#i Updating seirPackage documentation
#i Loading seirPackage
