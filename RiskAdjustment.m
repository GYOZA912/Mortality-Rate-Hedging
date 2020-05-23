function [SurvFixedOURA,SurvFixedFELRA,SurvFixedOU,SurvFixedFEL,PriceOURA,PriceFELRA,LambdaFEL,LambdaOU] = ...
    RiskAdjustment(qParamOU,qParamFEL,NbPeriod,RiskFactor,ActuFactorOU,ActuFactorFEL,q0)

LambdaMCOU = -0.0071;
LambdaMCFEL = -0.0071;
T = 1 : NbPeriod;

%%
Express1 = ((qParamOU(1,2)^2)/(2*qParamOU(1,1)^2))*T;
Express2 = ((qParamOU(1,2)^2)/(qParamOU(1,1)^3))*exp(qParamOU(1,1)*T);
Express3 = ((qParamOU(1,2)^2)/(4*qParamOU(1,1)^3))*exp(2*qParamOU(1,1)*T);
Express4 = (3*qParamOU(1,2)^2)/(4*qParamOU(1,1)^3);

alpha_OU = Express1 - Express2 + Express3 + Express4;
beta_OU = (1 - exp(qParamOU(1,1)*T)) / qParamOU(1,1);

b = -sqrt((qParamFEL(1,1)^2) + 2 * (qParamFEL(1,2)^2));
c = (b + qParamFEL(1,1)) / 2;
d = (b - qParamFEL(1,1)) / 2;
beta_FEL = (1 - exp(b * T)) ./ (c + (d * exp(b * T)));
    
SurvFixedOU(1,T) = exp(alpha_OU + beta_OU .* q0);
SurvFixedFEL(1,T) = exp(beta_FEL .* q0); 

PriceOU = sum(SurvFixedOU(1,T).*RiskFactor(T).*ActuFactorOU(T));
PriceFEL = sum(SurvFixedFEL(1,T).*RiskFactor(T).*ActuFactorFEL(T));

%%
BESurvOURA = @(p) exp((((qParamOU(1,2)^2)/(2*(qParamOU(1,1) + p)^2)) * T) - ...
    (((qParamOU(1,2)^2)/((qParamOU(1,1) + p)^3))*exp((qParamOU(1,1)+p) * T))...
    + (((qParamOU(1,2)^2)/(4*(qParamOU(1,1) + p)^3))*exp(2*(qParamOU(1,1) + p) * T))...
    + ((3*qParamOU(1,2)^2)/(4*(qParamOU(1,1) + p)^3))...
    + ((1-exp((qParamOU(1,1) + p)*T))/(qParamOU(1,1) + p)) * q0);
objective =  @(p) sum((PriceOU - sum(BESurvOURA(p).*ActuFactorOU(T))).^2)/size(T,2);
LambdaOU = fmincon(objective,LambdaMCOU);

SurvFixedOURA = BESurvOURA(LambdaOU);

BESurvFELRA = @(y) exp((1-exp((-sqrt(((qParamFEL(1,1)+y)^2)+2*(qParamFEL(1,2)^2))) .* T))./...
    ((((-sqrt(((qParamFEL(1,1)+y)^2)+2*(qParamFEL(1,2)^2))) + (qParamFEL(1,1)+y))/2)...
    + (((-sqrt(((qParamFEL(1,1)+y)^2)+2*(qParamFEL(1,2)^2))) - (qParamFEL(1,1)+y))/2) ...
    * exp((-sqrt(((qParamFEL(1,1)+y)^2)+2*(qParamFEL(1,2)^2))) * T)) * q0);
objectiveFEL =  @(y) sum((PriceFEL - sum(BESurvFELRA(y).*ActuFactorFEL(T))).^2)/size(T,2);
LambdaFEL = fmincon(objectiveFEL,LambdaMCFEL);

SurvFixedFELRA = BESurvFELRA(LambdaFEL);

%%
PriceOURA = sum(BESurvOURA(LambdaOU).*ActuFactorOU(T));
PriceFELRA = sum(BESurvFELRA(LambdaFEL).*ActuFactorFEL(T));



end
