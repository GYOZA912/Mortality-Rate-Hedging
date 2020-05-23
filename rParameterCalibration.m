function [ActuFactorOU,ActuFactorFEL,rParamOU,rParamFEL,r0] =rParameterCalibration(NbPeriod)

%%%================interest rate calbration===============%%%

    %%% simulate zero-coupon rate by finding a_hat theta_hat sigma_hat that minimize MSE
    %%%set original parameters for fmincon function
    %dataOriginal=xlsread('CADZC2018-12-21',1);%zero-coupon rate for 14 days of different years
    % for quarterly data

    load('interestrateUS.mat');
    TimeStep = 1;
    %r = dataOriginal(2,:);%if wanting to use annual ones,comment this line out and use following paragraph commented out instad
    r = A(:,2:end);
    r0 = A(:,1);
    %{
    for i = 1:30 %for annualy data
        r(1,i) = dataOriginal(2, 4*i);
    end
    %} 
    S_hist = exp(-cumsum(r,2)*TimeStep);
    lb = [0,0,0];
    ub = [3,1,1];
    p0 = [1.5,0.01,0.01];
    
    %%%simulation zero-coupon rate with formula S(t) = exp[alpha(t) + beta(t)*r0]
    n = 1: size(r,2); 
    X = 1: size(r,1); 
    t = n;
    T = 1 : NbPeriod;
    
    
    %% Vasicek model
        %%%pV(1) = k, pV(2) = gamma, pV(3) = sigma
        betaVasicek = @(pV) (exp(-pV(1)*t) - 1)/pV(1);
        alphaVasicek = @(pV) -(betaVasicek(pV)+t) * (pV(1)^2*pV(2)-pV(3)^2/2)/pV(1)^2 - pV(3)^2*betaVasicek(pV).^2/(4*pV(1));
        S_Vasicek = @(pV) exp(alphaVasicek(pV) + betaVasicek(pV).*r0);
        ObjectVasicek = @(pV) mean(sum((S_hist(X,t) - S_Vasicek(pV)).^2)); %MSE of simulated and real zero-coupon bond price as objectve function
    
        [rParamOU,MseVasicek] = fmincon(ObjectVasicek,p0,[],[],[],[],lb,ub);
        
        %simulate for years = NbPeriod
        betaVasicek_T = @(pV) (exp(-pV(1)*T) - 1)/pV(1);
        alphaVasicek_T = @(pV) -(betaVasicek_T(pV)+T) * (pV(1)^2*pV(2)-pV(3)^2/2)/pV(1)^2 - pV(3)^2*betaVasicek_T(pV).^2/(4*pV(1));
        S_Vasicek_T = @(pV) exp(alphaVasicek_T(pV) + betaVasicek_T(pV)*r0(1));
        
        
    %% CIR model
        %%%pC(1) = k, pC(2) = gamma, pC(3) = sigma
        b = @(pC) -sqrt(pC(1)^2+2*pC(3)^2);
        c = @(pC) (b(pC)- pC(1))/2;
        d = @(pC) (b(pC) + pC(1))/2;
        alphaCIR = @(pC) -(2*pC(1)*pC(2))/pC(3)^2 * log((c(pC)+d(pC)*exp(b(pC).*t))/b(pC)) + pC(1)*pC(2)/c(pC).*t;
        betaCIR = @(pC) (1-exp(b(pC).*t))/(c(pC)+d(pC)*exp(b(pC).*t));
        S_Cir = @(pC) exp(alphaCIR(pC) + betaCIR(pC)*r0);
        ObjectCir = @(pC) mean(sum((S_hist(X,t) - S_Cir(pC)).^2)); %MSE of simulated and real zero-coupon bond price as objectve function
    
        [rParamFEL,MseCir] = fmincon(ObjectCir,p0,[],[],[],[],lb,ub);
    
        %For simulate for years = NbPeriod
        alphaCIR_T = @(pC) -(2*pC(1)*pC(2))/pC(3)^2 * log((c(pC)+d(pC)*exp(b(pC).*T))/b(pC)) + pC(1)*pC(2)/c(pC).*T;
        betaCIR_T = @(pC) (1-exp(b(pC).*T))/(c(pC)+d(pC)*exp(b(pC).*T));
        S_Cir_T = @(pC) exp(alphaCIR_T(pC) + betaCIR_T(pC)*r0(1));
        
        
    %% get simulated zero-coupon bond price and rate 
    ActuFactorOU = S_Vasicek_T(rParamOU);
    ActuFactorFEL = S_Cir_T(rParamFEL);
    
    %% plot
    subplot(1,2,1)
    plot(66:100, Survival1966M(1,end-NbPeriod+1:end));
    hold on; 
    plot(66:100, SurvFixedOU);
    hold on
    plot(66:100, SurvFixedFEL);
    legend('历史记录','OU','Feller')
    xlael('年龄') 
    ylabel('生b存率') 
    title('生存率模拟')
    
    subplot(1,2,2)
    plot(1:29,S_hist(1,1:29));
    hold on;
    plot(1:29,ActuFactorOU(:,1:29));
    hold on;
    plot(1:29,ActuFactorFEL(:,1:29));
    xlabel('日历时间') 
    ylabel('零息债券价格')
    legend('历史数据','Vesicek','Cir')
    title('零息债券价格模拟')
end