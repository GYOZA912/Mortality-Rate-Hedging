
function [SurplusOU,SurplusFEL,SwapSurplusOU,SwapSurplusFEL,CapSurplusOU,CapSurplusFEL,...%KQDSurplusOU,KQDSurplusFEL,...
    Payoff_Swap_OU,Payoff_Swap_FEL,Payoff_Cap_OU,Payoff_Cap_FEL] ...
    = HedgingSim(NbSims,NumKeyAge,KeyAgeMin,KeyAgeMax,PriceOURA,PriceFELRA,qParamOU,qParamFEL,...
    ActuFactorOU,ActuFactorFEL,SurvFixedOU,SurvFixedFEL,SurvFixedOURA,SurvFixedFELRA,...
    rParamOU,rParamFEL,LambdaOU,LambdaFEL,Nb_Annuitants,NbPeriod,DureeContract,cohort,q0,r0)


%%Assets Total Value 
    AssetsOU = PriceOURA*Nb_Annuitants;%in total 4000 annuitants
    AssetsFEL = PriceFELRA*Nb_Annuitants;
%{
    qBE_OURA = (1+qParamOU(1)+LambdaOU).^(1:NbPeriod)*q0;
    qBE_FELRA = (1+qParamFEL(1)+LambdaFEL).^(1:NbPeriod)*q0;
    qBE_OU = (1+qParamOU(1)).^(1:NbPeriod)*q0;
    qBE_FEL = (1+qParamFEL(1)).^(1:NbPeriod)*q0;
    
    AgeInterval = (KeyAgeMax - KeyAgeMin)/ (NumKeyAge-1);
    [q_jOU,q_jFEL] = KQD(qBE_OU,qBE_FEL,NumKeyAge,AgeInterval,KeyAgeMax,KeyAgeMin,NbPeriod,cohort);
%}
    
for i = 1: NbSims
%% Monte Carlo for interest rate
    %%%get mean for interest rate from 100 simulations
    for j = 1:50
        %%%dt is omitted since it's 1
        
        %The first values of both model are controlled to be larger than 0
        SimRateOU(j,1) = 0;
        SimRateFEL(j,1) = 0;
        while SimRateOU(j,1)<=0
        
            SimRateOU(j,1) = r0(1) + rParamOU(1)*...
            (rParamOU(2)-r0(1)) + rParamOU(3)*randn;
        end    
        while SimRateFEL(j,1)<=0   
           SimRateFEL(j,1) = r0(1) + rParamFEL(1)*...
           (rParamFEL(2)-r0(1)) + sqrt(max(r0(1),0))*rParamFEL(3)*randn;
        end

        for Period = 2 : NbPeriod 
            SimRateOU(j,Period) = SimRateOU(j,Period-1) + rParamOU(1)*...
                (rParamOU(2)-SimRateOU(j,Period-1)) + rParamOU(3)*randn;            
            while SimRateOU(j,Period) <= 0            
                SimRateOU(j,Period) = SimRateOU(j,Period-1) + rParamOU(1,1)*...
                (rParamOU(1,2)-SimRateOU(j,Period-1)) + rParamOU(1,3)*randn;
            end
            

            SimRateFEL(j,Period) = SimRateFEL(j,Period-1) + rParamFEL(1)*...
                (rParamFEL(2)-SimRateFEL(j,Period-1)) + sqrt(max(SimRateFEL(j,Period-1),0))*rParamFEL(3)*randn;
            while SimRateFEL(1,Period) <= 0            
                SimRateFEL(1,Period) = SimRateFEL(1,Period-1) + rParamFEL(1,1)*...
                (rParamFEL(1,2)-SimRateFEL(1,Period-1)) + rParamFEL(1,3)*randn;
            end

        end
    end

    MeanSimRateOU(i,:) = mean(SimRateOU,1);
    MeanSimRateFEL(i,:) = mean(SimRateFEL,1);

    ActuFactorOUSim(i,:) = exp(-cumsum(MeanSimRateOU(i,:),2));
    ActuFactorFELSim(i,:) = exp(-cumsum(MeanSimRateFEL(i,:),2));


%% Monte for mortality rate

    for T = 1 : NbPeriod
    
        Alea = randn(1,NbPeriod); %%% Random normal independant variables


        if T == 1
                ForecastOU(1,T) = q0*(1+qParamOU(1))+qParamOU(2).*Alea(1,1); 
                ForecastFEL(1,T) = q0*(1+qParamFEL(1))+qParamFEL(2)*sqrt(max(q0,0)).*Alea(1,1); 
                ForecastOURA(1,T) = q0 * (1+qParamOU(1,1)+LambdaOU) + qParamOU(1,2)*Alea(1,1);
                ForecastFELRA(1,T) = q0 * (1+qParamFEL(1,1)+LambdaFEL) + qParamFEL(1,2)*sqrt(max(0,q0))*Alea(1,1);
        else
                ForecastOU(1,T) = ForecastOU(1,T-1)*(1+qParamOU(1))+qParamOU(2).*Alea(1,T); 
                ForecastFEL(1,T) = ForecastFEL(1,T-1)*(1+qParamFEL(1))+qParamFEL(2)*sqrt(max(ForecastFEL(1,T-1),0)).*Alea(1,T); 
                ForecastOURA(1,T) = ForecastOURA(1,T-1)* (1+qParamOU(1,1)+LambdaOU) + qParamOU(1,2)*Alea(1,T);
                ForecastFELRA(1,T) = ForecastFELRA(1,T-1)* (1+qParamFEL(1,1)+LambdaFEL) + qParamFEL(1,2)*sqrt(max(0,ForecastFEL(1,T-1)))*Alea(1,T);
        end 

    end


%% Simulate death time

    
    RandomExp = -log(1-rand(1,Nb_Annuitants)); %%% Variable aléatoire exponentielle
    
    ScaleOU = cumsum(ForecastOU,2); 
    ScaleFEL = cumsum(ForecastFEL,2); 
    %ScaleOU_KQD = cumsum(q_jOU,2); 
    %ScaleFEL_KQD = cumsum(q_jFEL,2);  
    
    for T = 1 : NbPeriod
       %{
        for j = 1:NumKeyAge 
        	DT_OU_KQD = RandomExp < ScaleOU_KQD(j,T);
            DT_FEL_KQD = RandomExp < ScaleFEL_KQD(j,T);
            LiabilitiesPerPeriodOU_KQD(j,T) = (Nb_Annuitants - sum(DT_OU_KQD,2))*ActuFactor(1,T); 
            LiabilitiesPerPeriodFEL_KQD(j,T) = (Nb_Annuitants - sum(DT_FEL_KQD,2))*ActuFactor(1,T);
            AnnuitantsLeftOU_KQD(j,T) = Nb_Annuitants-sum(DT_OU_KQD,2);
            AnnuitantsLeftFEL_KQD(j,T) = Nb_Annuitants-sum(DT_FEL_KQD,2);
        end
       %}
        
        DT_OU = RandomExp < ScaleOU(1,T); %%% : 1 x 4000 (Where True = 1 = Décès)
        DT_FEL = RandomExp < ScaleFEL(1,T);      
        
        LiabilitiesPerPeriodOU(1,T) = (Nb_Annuitants - sum(DT_OU,2))*ActuFactorOUSim(1,T); 
        LiabilitiesPerPeriodFEL(1,T) = (Nb_Annuitants - sum(DT_FEL,2))*ActuFactorFELSim(1,T); 
        
        RealizedSurvieOU(1,T) = 1-(sum(DT_OU,2)/Nb_Annuitants);
        RealizedSurvieFEL(1,T) = 1-(sum(DT_FEL,2)/Nb_Annuitants);
       
    
    end  


%% Calculate returns
    %%% unhedged
    TotalLiabilitiesOU = sum(LiabilitiesPerPeriodOU,2)';
    TotalLiabilitiesFEL = sum(LiabilitiesPerPeriodFEL,2)';
    SurplusOU(i,:) =(AssetsOU - TotalLiabilitiesOU)/Nb_Annuitants; %%Surplus per contract
    SurplusFEL(i,:) =(AssetsFEL - TotalLiabilitiesFEL)/Nb_Annuitants; 
    
    %{
    %%%KQD

    for j = 1:NumKeyAge
        Tj=AgeInterval*(j-1)+KeyAgeMin-cohort;
        
        TotalLiabilitiesOU_KQD(j) = sum(LiabilitiesPerPeriodOU_KQD(j,:),2)';
        TotalLiabilitiesFEL_KQD(j) = sum(LiabilitiesPerPeriodFEL_KQD(j,:),2)';
    
        KQD_Pj_OU(j) = (TotalLiabilitiesOU - TotalLiabilitiesOU_KQD(j));
        KQD_Pj_FEL(j) = (TotalLiabilitiesFEL - TotalLiabilitiesFEL_KQD(j)); 
        WeightOU(j) = KQD_Pj_OU(j)/(-(1+InterestRate)^(-Tj));
        WeightFEL(j) = KQD_Pj_FEL(j)/(-(1+InterestRate)^(-Tj));
        
        qRealizedOU(j) = 1-AnnuitantsLeftOU_KQD(j,Tj)/AnnuitantsLeftOU_KQD(j,Tj-1);
        qRealizedFEL(j) = 1-AnnuitantsLeftFEL_KQD(j,Tj)/AnnuitantsLeftFEL_KQD(j,Tj-1);
        
        HedgeOU_KQD(i,j) = WeightOU(j)*((2*qBE_OU(Tj)-qBE_OURA(Tj))-qRealizedOU(j))/(1+InterestRate)^(Tj);
        HedgeFEL_KQD(i,j) = WeightFEL(j)*((2*qBE_FEL(Tj)-qBE_FELRA(Tj))-qRealizedFEL(j))/(1+InterestRate)^(Tj);

        %HedgeOU_KQD(j) = WeightOU(j)*(qRealizedOU(j)-(2*qBE_OU(AgeInterval*j)-ForecastOURA(AgeInterval*j)))/(1+InterestRate)^(AgeInterval*j);
        %HedgeFEL_KQD(j) = WeightFEL(j)*(qRealizedFEL(j)-(2*qBE_FEL(AgeInterval*j)-ForecastFELRA(AgeInterval*j)))/(1+InterestRate)^(AgeInterval*j);
        
        %HedgeOU_KQD(j) = WeightOU(j)*(qRealizedOU(j)-qBE_OU(AgeInterval*j))/(1+InterestRate)^(AgeInterval*j);
        %HedgeFEL_KQD(j) = WeightFEL(j)*(qRealizedFEL(j)-SurvFixedFEL(AgeInterval*j))/(1+InterestRate)^(AgeInterval*j);
    end
    %}
    %%% hedged
    ActuFactorOUkDuree = ActuFactorOUSim(:,1:DureeContract)';
    ActuFactorFELDuree = ActuFactorFELSim(:,1:DureeContract)';
    
    Payoff_Swap_OU(i,:) = sum(ActuFactorOUkDuree(:,1).*(RealizedSurvieOU(1,1:DureeContract) - SurvFixedOURA(1,1:DureeContract))');
    Payoff_Swap_FEL(i,:)= sum(ActuFactorFELDuree(:,1).*(RealizedSurvieFEL(1,1:DureeContract) - SurvFixedFELRA(1,1:DureeContract))');     

    Payoff_Cap_OU(i,:) = sum(ActuFactorOUkDuree(:,1).*(max(RealizedSurvieOU(1,1:DureeContract) - SurvFixedOU(1,1:DureeContract),0))'); 
    Payoff_Cap_FEL(i,:) = sum(ActuFactorFELDuree(:,1).*(max(RealizedSurvieFEL(1,1:DureeContract) - SurvFixedFEL(1,1:DureeContract),0))'); 



end

%{
KQDSurplusOU = SurplusOU + sum(HedgeOU_KQD,2)/Nb_Annuitants;
KQDSurplusFEL = SurplusFEL + sum(HedgeFEL_KQD,2)/Nb_Annuitants;  
%}

SwapSurplusOU = SurplusOU + Payoff_Swap_OU;
SwapSurplusFEL = SurplusFEL + Payoff_Swap_FEL;
    
[C_Cap_OU,C_Cap_FEL] = CapPricing(SurvFixedOURA,SurvFixedFELRA,SurvFixedOU,SurvFixedFEL,Nb_Annuitants,DureeContract,...
    NbSims,LambdaOU,LambdaFEL,NbPeriod,qParamOU,qParamFEL,q0,ActuFactorOU,ActuFactorFEL);
CapSurplusOU = SurplusOU + Payoff_Cap_OU - C_Cap_OU;
CapSurplusFEL = SurplusFEL + Payoff_Cap_FEL - C_Cap_FEL;

end