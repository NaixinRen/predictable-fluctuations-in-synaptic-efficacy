function b1 = random_parameter_s2(y,XX,t,b_s1,b_glm,distance,v,eta,tau0,stage)

% t = linspace(-25,25,102);
% t = t+mean(diff(t))/2;
% t = t(1:101);


options=[];
options.method = 'cg';
options.MaxIter = 50;
options.Display = 'off';
f=Inf;


dt_st = min(distance*v(1)+v(2),max(t)-1);
synt = t;
synt(t<dt_st)=dt_st;
syn = stage*(synt-dt_st)/tau0.*exp(1-(synt-dt_st)/tau0);

b_est= glmfit([XX',syn(1:length(y))'],y','poisson','constant','off');
w_st = abs(b_est(end))+(b_est(end)==0)*.001;

switch stage
    
    case 1
        
        for rr = 1:4
            %             if mod(rr,2) == 0
            %                 b0 = [b_est(1:end-1);log(w_st); log(dt_st); log(tau0)];
            %             else
            %                 b0 = b_s1';
            %             end
            
            switch mod(rr,4)
                case 0
                    b0 = [b_glm';log(.001);log(dt_st); log(tau0)];
                case 1
                    b0 = [b_glm(1); b_glm(2:end)'+randn(size(XX,1)-1,1)/5;log(w_st);log(dt_st); log(tau0)];
                case 2
                    b0 = [log(nanmean(y)); b_glm(2:end)'*0;log(w_st); log(dt_st); log(tau0)];
                case 3
                    b0 = b_s1';
            end
            
            [brr,frr] = minFunc(@loss_excalpha,b0,options,XX',y',t',v,distance,eta,tau0);
            if frr<f
                b1=brr;
                f=frr;
            end
        end
        
    case -1
        
        for rr = 1:4
            %             if mod(rr,2) == 0
            %                 b0 = [b_est(1:end-1);log(w_st); log(dt_st); log(tau0)];
            %             else
            %                 b0 = b_s1';
            %             end
            switch mod(rr,4)
                case 0
                    b0 = [b_glm';log(.001);log(dt_st); log(tau0)];
                case 1
                    b0 = [b_glm(1); b_glm(2:end)'+randn(size(XX,1)-1,1)/5;log(w_st);log(dt_st); log(tau0)];
                case 2
                    b0 = [log(nanmean(y)); b_glm(2:end)'*0;log(w_st); log(dt_st); log(tau0)];
                case 3
                    b0 = b_s1';
            end
            [brr,frr] = minFunc(@loss_inhalpha,b0,options,XX',y',t',v,distance,eta,tau0);
            if frr<f
                b1=brr;
                f=frr;
            end
        end
        
end

