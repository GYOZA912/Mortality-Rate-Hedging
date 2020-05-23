function [HedgeSingle65Mean,HedgeSingle65Std,HedgeSingle65Skew,HedgeSingle65Kutosis,...
    HedgeSingle65VaR,HedgeSingle65ES] = AllTests(SurplusOU,SurplusFEL,SwapSurplusOU,...
    SwapSurplusFEL,CapSurplusOU,CapSurplusFEL)%,KQDSurplusOU,KQDSurplusFEL)
%KQDSurplusOU_new=KQDSurplusOU(~isnan(KQDSurplusOU));
%KQDSurplusFEL_new=KQDSurplusFEL(~isnan(KQDSurplusFEL));

%%%Mean
HedgeSingle65Mean(1,1) = mean(SurplusOU);
HedgeSingle65Mean(2,1) = mean(SurplusFEL);
HedgeSingle65Mean(3,1) = mean(SwapSurplusOU);
HedgeSingle65Mean(4,1) = mean(SwapSurplusFEL);
HedgeSingle65Mean(5,1) = mean(CapSurplusOU);
HedgeSingle65Mean(6,1) = mean(CapSurplusFEL);
%HedgeSingle65Mean(7,1) = mean(KQDSurplusOU_new);
%HedgeSingle65Mean(8,1) = mean(KQDSurplusFEL_new);
%%%std
HedgeSingle65Std(1,1) = std(SurplusOU);
HedgeSingle65Std(2,1) = std(SurplusFEL);
HedgeSingle65Std(3,1) = std(SwapSurplusOU);
HedgeSingle65Std(4,1) = std(SwapSurplusFEL);
HedgeSingle65Std(5,1) = std(CapSurplusOU);
HedgeSingle65Std(6,1) = std(CapSurplusFEL);
%HedgeSingle65Std(7,1) = std(KQDSurplusOU_new);
%HedgeSingle65Std(8,1) = std(KQDSurplusFEL_new);
%%%skewness
HedgeSingle65Skew(1,1) = skewness(SurplusOU);
HedgeSingle65Skew(2,1) = skewness(SurplusFEL);
HedgeSingle65Skew(3,1) = skewness(SwapSurplusOU);
HedgeSingle65Skew(4,1) = skewness(SwapSurplusFEL);
HedgeSingle65Skew(5,1) = skewness(CapSurplusOU);
HedgeSingle65Skew(6,1) = skewness(CapSurplusFEL);
%HedgeSingle65Skew(7,1) = skewness(KQDSurplusOU_new);
%HedgeSingle65Skew(8,1) = skewness(KQDSurplusFEL_new);
%%%kurtosis
HedgeSingle65Kutosis(1,1) = kurtosis(SurplusOU);
HedgeSingle65Kutosis(2,1) = kurtosis(SurplusFEL);
HedgeSingle65Kutosis(3,1) = kurtosis(SwapSurplusOU);
HedgeSingle65Kutosis(4,1) = kurtosis(SwapSurplusFEL);
HedgeSingle65Kutosis(5,1) = kurtosis(CapSurplusOU);
HedgeSingle65Kutosis(6,1) = kurtosis(CapSurplusFEL);
%HedgeSingle65Kutosis(7,1) = kurtosis(KQDSurplusOU_new);
%HedgeSingle65Kutosis(8,1) = kurtosis(KQDSurplusFEL_new);

%%%Var & BE
VaRLevel = 0.99;

[VaR,ES] = VaRES(SurplusOU,VaRLevel);
HedgeSingle65VaR(1,1) = VaR;
HedgeSingle65ES(1,1) = ES;
[VaR,ES] = VaRES(SurplusFEL,VaRLevel);
HedgeSingle65VaR(2,1) = VaR;
HedgeSingle65ES(2,1) = ES;
[VaR,ES] = VaRES(SwapSurplusOU,VaRLevel);
HedgeSingle65VaR(3,1) = VaR;
HedgeSingle65ES(3,1) = ES;
[VaR,ES] = VaRES(SwapSurplusFEL,VaRLevel);
HedgeSingle65VaR(4,1) = VaR;
HedgeSingle65ES(4,1) = ES;
[VaR,ES] = VaRES(CapSurplusOU,VaRLevel);
HedgeSingle65VaR(5,1) = VaR;
HedgeSingle65ES(5,1) = ES;
[VaR,ES] = VaRES(CapSurplusFEL,VaRLevel);
HedgeSingle65VaR(6,1) = VaR;
HedgeSingle65ES(6,1) = ES;
%[VaR,ES] = VaRES(KQDSurplusOU_new,VaRLevel);
%HedgeSingle65VaR(7,1) = VaR;
%HedgeSingle65ES(7,1) = ES;
%[VaR,ES] = VaRES(KQDSurplusFEL_new,VaRLevel);
%HedgeSingle65VaR(8,1) = VaR;
%HedgeSingle65ES(8,1) = ES;

end