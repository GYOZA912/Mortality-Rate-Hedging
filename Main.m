clear

%% initialization
load('outputCHNfemales.mat')
LengthDevelopment = 20;
cohort = 70; 
MaxAge = 90;
NbPeriod = 20;%MaxAge - cohort;
Nb_Annuitants = 4000; %%% Number of rentiers, by cohort, contained in the simulated portfolio
NbSims = 100000;%for MC simulation in hedging and premium
RiskPrem = 0.002;
%InterestRate = 0.04;
DureeContract = 15;
RiskFactor = exp(RiskPrem*(1:NbPeriod));
%ActuFactor = exp(-InterestRate*(1:NbPeriod)); %%% Zero-coupon prices : [1xNbPeriod]
NbGroup = 1; %For single-cohort simulation
NumKeyAge = 6;
KeyAgeMin = 70;
KeyAgeMax = 95;

%% Create Survival Matrix
%%% read history data of corhort of 65 from 1966-2016

output = -log(1-(output./(1+0.5.*output)));
Data2016 = output(end - LengthDevelopment: end,end - LengthDevelopment : end); %%% 1966 - 2016

for x = 1:size(Data2016,1)
    for t = 1:size(Data2016,2)

        if t == 1 || x == 1
                Survival1966(x,t) = 1 - Data2016(x,t);
        else
                Survival1966(x,t) = Survival1966(x-1,t-1)*(1 - Data2016(x,t));
        end

    end
end

for x = 1:size(Data2016,1)
    for t = 1:size(Data2016,2)

        if x + t - 1 > size(Data2016,1)
                break
        elseif t == 1
                Data1966_History(x,t) = Data2016(x,t);%%Record mortality rate
                Survival1966M(x,t) = Survival1966(x,t);%%Record mortality rate
        else
                Data1966_History(x,t) = Data2016(x+t-1,t);
                Survival1966M(x,t) = Survival1966(x+t-1,t);
        end

    end
end

q0 = Data1966_History(1,end-(MaxAge-cohort));


%% calibrations
%%%Interest rate calibration
[ActuFactorOU,ActuFactorFEL,rParamOU,rParamFEL,r0] =rParameterCalibration(NbPeriod);

%%%mortality rate calibration
[qParamOU,qParamFEL]= qParameterCalibration(Data1966_History,Survival1966M);

[SurvFixedOURA,SurvFixedFELRA,SurvFixedOU,SurvFixedFEL,PriceOURA,PriceFELRA,LambdaFEL,LambdaOU] = ...
    RiskAdjustment(qParamOU,qParamFEL,NbPeriod,RiskFactor,ActuFactorOU,ActuFactorFEL,q0);
    
%% Monte carlo for portfolio simulation
[SurplusOU,SurplusFEL,SwapSurplusOU,SwapSurplusFEL,CapSurplusOU,CapSurplusFEL,...%KQDSurplusOU,KQDSurplusFEL,...
    Payoff_Swap_OU,Payoff_Swap_FEL,Payoff_Cap_OU,Payoff_Cap_FEL] ...
    = HedgingSim(NbSims,NumKeyAge,KeyAgeMin,KeyAgeMax,PriceOURA,PriceFELRA,qParamOU,qParamFEL,...
	ActuFactorOU,ActuFactorFEL,SurvFixedOU,SurvFixedFEL,SurvFixedOURA,SurvFixedFELRA,...
    rParamOU,rParamFEL,LambdaFEL,LambdaOU,Nb_Annuitants,NbPeriod,DureeContract,cohort,q0,r0);


%% tests
[HedgeSingle65Mean,HedgeSingle65Std,HedgeSingle65Skew,HedgeSingle65Kutosis,...
    HedgeSingle65VaR,HedgeSingle65ES] = AllTests(SurplusOU,SurplusFEL,SwapSurplusOU,...
    SwapSurplusFEL,CapSurplusOU,CapSurplusFEL);%,KQDSurplusOU,KQDSurplusFEL);
WayOfHedge = ["UN_OU";"UN_FEL";"SWAP_OU";"SWAP_FEL";"CAP_OU";"CAP_FEL"];%"KQD_OU";"KQD_FEL"];

TestsTable = table(WayOfHedge,HedgeSingle65Mean,HedgeSingle65Std,HedgeSingle65Skew,...
    HedgeSingle65Kutosis,HedgeSingle65VaR,HedgeSingle65ES);
