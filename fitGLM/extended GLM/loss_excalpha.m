function [f,df,lam,syn] = loss_excalpha(b,ACCG,X,y,t,v,dis,eta,tau0)

p = b((size(X,2)+1):end);
p=exp(p);
b = b(1:size(X,2));

w = p(1);
deltat = p(2); 
tau = p(3);
eta_w = eta(1);
eta_dt = eta(2);
eta_tau = eta(3);

t(t<deltat)=deltat;
syn = w*(t-deltat)/tau.*exp(1-(t-deltat)/tau);
syn = conv(syn,ACCG/max(ACCG),'same');
mm = max(abs(syn))+ (max(abs(syn))==0);
syn = syn./mm*abs(w);

lam = exp(X*b + syn);
f = -nansum((y.*log(lam+(lam==0))-lam));
fw = eta_w*(w-0)^2;
fdt = eta_dt*((deltat)-dis*v(1)-v(2))^2; 
ftau = eta_tau*(tau-.8)^2;
f = f+fw+fdt+ftau;

dl = (lam - y);
dl(~isfinite(dl)) = 0;
db = X'*dl;
dw = syn'*dl/w+eta_w*2*(w-0); 
tmp = (syn/tau - syn./(t-deltat)).*dl;
tmp(isinf(tmp)) = nan;
ddeltat = nansum(conv(tmp(~isnan(tmp)),ACCG/max(ACCG),'same'))+eta_dt*2*((deltat)-dis*v(1)-v(2));
tmp = ((-syn/tau+syn.*(t-deltat)/tau.^2)).*dl;
dtau = nansum(conv(~isnan(tmp), ACCG/max(ACCG),'same'))+eta_tau*2*(tau-tau0);
% ddeltat = nansum(tmp)+eta_dt*2*((deltat)-dis*v(1)-v(2));
% dtau = nansum(((-syn/tau+syn.*(t-deltat)/tau.^2)).*dl)+eta_tau*2*(tau-tau0);
df =  [db;dw*w;ddeltat*deltat;dtau*tau]; 
syn = syn/abs(w);