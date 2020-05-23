function [C_Cap_OU,C_Cap_FEL] = CapPricing(SurvOUFixedRA,SurvFELFixedRA,SurvFixedOU,SurvFixedFEL,Nb_Annuitants,DureeContract,...
    NbSims,LambdaOU,LambdaFEL,NbPeriod,qParamOU,qParamFEL,q0,ActuFactorOU,ActuFactorFEL)

% Cette fontion utilise la moyenne des scénarios générés de taux de mortalité ajustés pour le risque afin de
% calculer la valeur moyenne des payoff actualisés d'un Cap ayant les strikes définis par Fung & al. (2015) 
%%%===================================%%%
%%%             INPUT
%%%===================================%%%

%%% Survival65RA : Moyenne des scénarios de taux de survies (Best estimate) de la cohorte dont on souhaite obtenir le prix du Cap
%%% Survival65 : Stikes du cap

%%%===================================%%%
%%%             OUTPUT
%%%===================================%%%

%%% C_Cap_Vas et C_Cap_Cir : Le prix d'une option

%%%===================================%%%
%[SurvOURA,SurvFELRA] = Simpath(NbSims,LambdaOU,LambdaFEL,NbPeriod,qParamOU,qParamFEL,q0);

 ActuFactorOUDuration = ActuFactorOU(:,1:DureeContract);
 ActuFactorFELDuration = ActuFactorFEL(:,1:DureeContract);

for Simul = 1:Nb_Annuitants
    
     RandomUni = -log(1-rand(1,Nb_Annuitants)); %%% Variable aléatoire exponentielle : 1 x 4000
     ScaleOU = -log(SurvOUFixedRA(1,1:DureeContract));%%%%%%here used fixed instead
     ScaleFEL = -log(SurvFELFixedRA(1,1:DureeContract));
 
    for Period = 1:DureeContract
        
        DTRA_OU = RandomUni < ScaleOU(1,Period);%%% True = 1 = Décès : 1 x 4000
        DTRA_FEL = RandomUni < ScaleFEL(1,Period);
        
        Obs_SurvieRAOU(Simul,Period) = 1-(sum(DTRA_OU,2)/Nb_Annuitants); 
        Obs_SurvieRAFEL(Simul,Period) = 1-(sum(DTRA_FEL,2)/Nb_Annuitants);  
    end
    
    D_Cap_OU(Simul,1:DureeContract) = max(ActuFactorOUDuration(1,:).*(Obs_SurvieRAOU(Simul,:) - SurvFixedOU(1,1:DureeContract)),0); 
    D_Cap_FEL(Simul,1:DureeContract) = max(ActuFactorFELDuration(1,:).*(Obs_SurvieRAFEL(Simul,:) - SurvFixedFEL(1,1:DureeContract)),0); 
    
end

C_Cap_OU = sum(mean(D_Cap_OU,1));
C_Cap_FEL = sum(mean(D_Cap_FEL,1));

end

