function[qParamOU,qParamFEL]= qParameterCalibration(Data1966_History,Survival1966M)

%% initialization
qParamOU = zeros(1,2);
qParamFEL = zeros(1,2);
MSEOU = zeros(1,1);
MSEFEL = zeros(1,1);

lb = [0,0];
ub = [1,1];
p0 = [0.1094,0.0007]; %%% Estimation based on male UK data born in 1945 (Luciano et Vigna (2008))
y0 = [0.11,0.007];

T = 1 : size(Data1966_History,2)-1;
X = 1 : size(Data1966_History,1)-1;

%% get optimal parameters for fixed survival rate

    %%%OU
    %%%p(1) for a, p(2) for sigma
    SurvivalFixedOU = @(p) max(exp(((p(2)^2)/(2*p(1)^2))*T-((p(2)^2)/(p(1)^3))*exp(p(1)*T)...
             +((p(2)^2)/(4*p(1)^3))*exp(2*p(1)*T)+((3*(p(2)^2))/(4*p(1)^3))...%%for alpha(t)
             +((1-exp(p(1)*T))/p(1)).*Data1966_History(X,1)),0);
    objectiveOU = @(p) mean(sum(sum((SurvivalFixedOU(p) - Survival1966M(X,T+1)).^2,2),1));     
    [qParamOU(1,:),MSEOU] = fmincon(objectiveOU,p0,[],[],[],[],lb,ub);
    
    %%%FELLER
    SurvivalFixedFEL = @(y) max(exp(((1-exp((-sqrt((y(1)^2)+(2*y(2)^2))).*T))./...
         ((((-sqrt((y(1)^2)+(2*y(2)^2)))+y(1))/2)+ (((-sqrt((y(1)^2)+(2*y(2)^2)))-y(1))/2)...
         *exp((-sqrt((y(1)^2)+(2*y(2)^2)))*T))).*Data1966_History(X,1)),0);
    objectiveFEL = @(y) mean(sum(sum((SurvivalFixedFEL(y) - Survival1966M(X,T+1)).^2,2),1));
    [qParamFEL(1,:),MSEFEL] = fmincon(objectiveFEL,y0,[],[],[],[],lb,ub);

    %%%get fixed survival rate for two models
    SurvFixedOU(X,T)= SurvivalFixedOU(qParamOU(1,:)');
    SurvFixedFEL(X,T)= SurvivalFixedFEL(qParamFEL(1,:)');
    


end