clear

%% 
load('CADfemaleslissage.mat')
cohort = 65; 
MaxAge = 100;
NbPeriod = MaxAge - cohort;
NbSims = 1000;%for MC simulation in hedging and premium

%%% read history data of corhort of 65 from 1982-2016
output = -log(1-(output./(1+0.5.*output)));
q0 = output(size(output,1) - NbPeriod ,size(output,2) - NbPeriod);
for i = size(output,1) - NbPeriod + 1 : size(output,1) %for row in original data
        Data65_History(1,i+NbPeriod-size(output,1)) = output(i,i+size(output,2)-size(output,1)); %for histoy data
end
Survival65 = exp(-cumsum(Data65_History,2));

qParamOU = zeros(size(Data65_History,1)-1,2);
qParamFEL = zeros(size(Data65_History,1)-1,2);
MSEOU = zeros(1,size(Data65_History,1)-1);
MSEFEL = zeros(1,size(Data65_History,1)-1);
SurvFixedOU = zeros(size(Data65_History,1),size(Data65_History,2));
SurvFixedFEL = zeros(size(Data65_History,1),size(Data65_History,2));

lb = [0,0.0001];
ub = [1,1];
p0 = [0.1094,0.0007]; %%% Estimation base on male UK data born in 1945 (Luciano et Vigna (2008))
y0 = [0.11,0.007];

T = 1 : 35;

%% get optimal parameters for fixed survival rate
for x = 1 : 1 %single cohort just 1        
    %%%OU
    %%%p(1) for a, p(2) for sigma, q0 for lambda(0)
    SurvivalFixedOU = @(p) max(exp(((p(2)^2)/(2*p(1)^2))*T-((p(2)^2)/(p(1)^3))*exp(p(1)*T)...
             +((p(2)^2)/(4*p(1)^3))*exp(2*p(1)*T)+((3*(p(2)^2))/(4*p(1)^3))...%%for alpha(t) part
             +((1-exp(p(1)*T))/p(1))*q0),0);         
    objectiveOU = @(p) mean((SurvivalFixedOU(p) - Survival65(x,T)).^2);     
    [qParamOU(x,:),MSEOU(x)] = fmincon(objectiveOU,p0,[],[],[],[],lb,ub);
    qProblemOUGlobal = createOptimProblem('fmincon','objective',objectiveOU,'x0',p0,'lb',lb,...
            'ub',ub,'options',optimset('Algorithm','SQP','Disp','none'));
    gs = GlobalSearch;
    qParamOUGlobal = run(gs,qProblemOUGlobal);
    
    %%%FELLER
    SurvivalFixedFEL = @(y) max(exp(((1-exp((-sqrt((y(1)^2)+(2*y(2)^2))).*T))./...
         ((((-sqrt((y(1)^2)+(2*y(2)^2)))+y(1))/2)+ (((-sqrt((y(1)^2)+(2*y(2)^2)))-y(1))/2)...
         *exp((-sqrt((y(1)^2)+(2*y(2)^2))).*T)))*q0),0);
    objectiveFEL = @(y) mean((SurvivalFixedFEL(y) - Survival65(x,T)).^2);
    [qParamFEL(x,:),MSEFEL(x)] = fmincon(objectiveFEL,y0,[],[],[],[],lb,ub);
    qproblemFELGlobal = createOptimProblem('fmincon','objective',objectiveFEL,'x0',y0,'lb',lb,...
                'ub',ub,'options',optimset('Algorithm','SQP','Disp','none'));
    qParamFELGlobal = run(gs,qproblemFELGlobal);
end
%%
for i = 1 : 10
    for Period = 2 : 35
    Alea = randn(1,35); %%% Random normal independant variables
    Forecast65FEL(i,1) = 0;
    Forecast65FEL_Global(i,1) = 0;

    Forecast65OU(i,1) = q0*(1+qParamOU(1))+qParamOU(2).*Alea(1,1); 

while Forecast65FEL(i,1)<=0
    Forecast65FEL(i,1) = q0*(1+qParamFEL(1))+qParamFEL(2)*sqrt(max(q0,0)).*Alea(1,1);
end

     Forecast65OU_Global(i,1) = q0*(1+qParamOUGlobal(1))+qParamOUGlobal(2).*Alea(1,1); 

while Forecast65FEL_Global(i,1)<=0
    Forecast65FEL_Global(i,1)= q0*(1+qParamFELGlobal(1))+qParamFELGlobal(2)*sqrt(max(q0,0)).*Alea(1,1);     
end
                Forecast65OU(i,Period) = q0*(1+qParamOU(1))+qParamOU(2).*Alea(1,1); 
                Forecast65FEL(i,Period) = q0*(1+qParamFEL(1))+qParamFEL(2)*sqrt(max(q0,0)).*Alea(1,1);
                Forecast65OU_Global(i,Period) = q0*(1+qParamOUGlobal(1))+qParamOUGlobal(2).*Alea(1,1); 
                Forecast65FEL_Global(i,Period) = q0*(1+qParamFELGlobal(1))+qParamFELGlobal(2)*sqrt(max(q0,0)).*Alea(1,1);     

                Forecast65OU(i,Period) = Forecast65OU(i,Period-1)*(1+qParamOU(1))+qParamOU(2).*Alea(1,Period); 
                Forecast65FEL(i,Period) = Forecast65FEL(i,Period-1)*(1+qParamFEL(1))+qParamFEL(2)*sqrt(max(Forecast65FEL(i,Period-1),0)).*Alea(1,Period); 
                Forecast65OU_Global(i,Period) = Forecast65OU_Global(i,Period-1)*(1+qParamOUGlobal(1))+qParamOUGlobal(2).*Alea(1,Period); 
                Forecast65FEL_Global(i,Period) = Forecast65FEL_Global(i,Period-1)*(1+qParamFELGlobal(1))+qParamFELGlobal(2)*sqrt(max(Forecast65FEL_Global(i,Period-1),0)).*Alea(1,Period);        
    end
end
MeanForecast65OU = mean(Forecast65OU,1);
MeanForecast65FEL = mean(Forecast65FEL,1);
MeanForecast65OU_Global = mean(Forecast65OU_Global,1);
MeanForecast65FEL_Global = mean(Forecast65FEL_Global,1);

S_OU = cat(2,Data65_History',MeanForecast65OU',MeanForecast65FEL'...
    ,MeanForecast65OU_Global',MeanForecast65FEL_Global');%merging three prices
figure(1)
plot(1:35,S_OU)
title('Real and Simulated Zero-coupon bond Price')
xlabel('Year')
ylabel('Rate')
legend({'Real q','OU Simulation','FEL Simulation','OU Simulation Gloabl','FEL Simulation Global'},...
    'location','Northeast')
