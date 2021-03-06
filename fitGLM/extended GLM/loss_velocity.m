function [f,df,err] = loss_velocity(V,dis,dt,weight,eta_delay)
v = exp(V(1));
delay = exp(V(2));
dt_est = dis*v+delay;
err = dt-dt_est;
f = nansum((dt-dt_est).^2.*weight)/sum(weight);

fdelay = eta_delay*(delay-0)^2;
f = f+fdelay;
dv = nansum(2*(dt-dt_est).*weight.*(-dis))/sum(weight);
ddelay = nansum(2*(dt-dt_est).*weight.*(-1))/sum(weight)+2*eta_delay*(delay-0);
df = [dv*v;ddelay*delay];