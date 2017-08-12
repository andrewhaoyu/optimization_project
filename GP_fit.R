library(GPfit)
n = 5; d = 1;
computer_simulator <- function(x){
  x = 2*x+0.5;
  y = sin(10*pi*x)/(2*x) + (x-1)^4;
  return(y)
}
set.seed(3);
library(lhs);
x = maximinLHS(n,d);
y = computer_simulator(x);
GPmodel = GP_fit(x,y);
print(GPmodel)
