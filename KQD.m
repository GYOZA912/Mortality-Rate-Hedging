function [q_jOU,q_jFEL] = KQD(qBE_OU,qBE_FEL,NumKeyAge,AgeInterval,KeyAgeMax,KeyAgeMin,NbPeriod,cohort)
KeyAge = KeyAgeMin:AgeInterval:KeyAgeMax;

%% hh
shift = zeros(NumKeyAge,NbPeriod);
q_jOU = zeros(NumKeyAge,NbPeriod);
q_jFEL = zeros(NumKeyAge,NbPeriod);
for j = 1:NumKeyAge

    if j == 1
        for T = 1:NbPeriod
            if cohort + T <= KeyAge(j)
                shift(j,T) = 0.001;
            else
                shift(j,T) = 0;
            end
        end
        
    elseif j ==NumKeyAge
        for T = 1:NbPeriod
            if cohort + T >= KeyAge(j)
                shift(j,T) = 0.001;
            else
                shift(j,T) = 0;
            end
        end
        
    else
        for T = 1:NbPeriod
            if cohort + T <= KeyAge(j-1) || cohort + T >= KeyAge(j+1)
                shift(j,T) = 0;
            elseif cohort + T <= KeyAge(j)
                shift(j,T) = 0.001*(cohort + T- KeyAge(j-1))/(KeyAge(j)-KeyAge(j-1));
            else
                shift(j,T) = 0.001*(KeyAge(j+1)-(cohort + T))/(KeyAge(j+1)-KeyAge(j));
            end
        end
    end
    
    q_jOU(j,:) = qBE_OU + shift(j,:);
    q_jFEL(j,:) = qBE_FEL + shift(j,:);
end

end