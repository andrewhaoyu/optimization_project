library(rstan);
stan.model.B <- "
data {

int<lower=0> N[12];
}

transformed data {

vector[15] ones;

for(i in 1:15){
ones[i] = 1;
}
}

parameters {

simplex[15] p;
real<lower=0, upper=1> q;

}

transformed parameters {

simplex[12] op;



op[1] = (1-q)*p[11];
op[2] = (1-q)*p[12];
op[3] = (1-q)*p[13];
op[4] = (1-q)*p[14];
op[5] = (1-q)*p[15];
op[6]= (1-q)*(p[1]+p[2]+p[3]+p[4]+p[5]+p[6]+p[7]+p[8]+p[9]+p[10]);
op[7] = q*(p[11]+p[6]);
op[8] = q*(p[12]+p[7]);
op[9] = q*(p[13]+p[8]);
op[10] = q*(p[14]+p[9]);
op[11] = q*(p[15]+p[10]);
op[12] = q*(p[1]+p[2]+p[3]+p[4]+p[5]);

}

model {
p  ~ dirichlet(ones);
q  ~ uniform(0,1);
N  ~ multinomial(op);
}

generated quantities {
real<lower=0, upper=5> tp[7];
real bounds[2];


tp[1] = (p[11]*1+p[12]*2+p[13]*3+p[14]*4+p[15]*5)/(p[11] + p[12]+p[13]+p[14]+p[15]);
tp[2] = p[11] + p[12] + p[13] + p[14] + p[15];
tp[3] = p[1]+p[2]+p[3]+p[4]+p[5] ;
tp[4] = p[6]+p[7]+p[8]+p[9]+p[10];
tp[5] = ((p[6]+p[11])+2*(p[7]+p[12])+3*(p[8]+p[13])+4*(p[9]+p[14])+5*(p[10]+p[15]))/(p[5]+p[6]+p[7]+p[8]+p[9]+p[10]+p[11]+p[12]+p[13]+p[14]+p[15]);
tp[6] = (p[6]+2*p[7]+3*p[8]+4*p[9]+5*p[10])/(p[6]+p[7]+p[8]+p[9]+p[10]);
tp[7] = 1*(p[1]+p[5]+p[11])+2*(p[2]+p[6]+p[12])+3*(p[3]+p[8]+p[13])+4*(p[4]+p[9]+p[13])+5*(p[5]+p[10]+p[15]);



bounds[1] = tp[5] * (tp[2] + tp[4]) + tp[3];
bounds[2] = tp[5] * (tp[2] + tp[4]) + 5*tp[3];
}

"

sm.bd <- stan_model(model_code =stan.model.B)
save("sm.bd",file="sm.bd.rdata")


stan.model.N <- "
data {

int<lower=0> N[6];
}

transformed data {

vector[10] ones;

for(i in 1:10){
ones[i] = 1;
}
}

parameters {

simplex[10] p;

}

transformed parameters {

simplex[6] op;



op[1] = p[1];
op[2] = p[2];
op[3] = p[3];
op[4] = p[4];
op[5] = p[5];
op[6]= p[6]+p[7]+p[8]+p[9]+p[10];

}

model {
p  ~ dirichlet(ones);
N  ~ multinomial(op);
}

generated quantities {

real bounds[2];
bounds[1] = p[1]+2*p[2]+3*p[3]+4*p[4]+5*p[5]+(p[6]+p[7]+
p[8]+p[9]+p[10]);
bounds[2] = p[1]+2*p[2]+3*p[3]+4*p[4]+5*p[5]+5*(p[6]+p[7]+
p[8]+p[9]+p[10]);
}

"


sm.na <- stan_model(model_code =stan.model.N)
save("sm.na",file="sm.na.rdata")



