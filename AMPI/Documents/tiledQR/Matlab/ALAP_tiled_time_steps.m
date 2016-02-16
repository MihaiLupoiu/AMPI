function [ Malap ] = ALAP_tiled_time_steps( Mcoarse )
%ALAP_TILED_TIME_STEPS Summary of this function goes here
%   Detailed explanation goes here

Malap = Mcoarse;
[rows, cols] = size(Malap);
ndiag = min(cols,rows);

for j=1:ndiag,
    if(j<rows)
        Mcoarse(j,j) = Mcoarse(j+1,j);
    end
    
    for i=j:rows,
        Malap(i,j) =  6*Mcoarse(2,1) + 22*(j-1) - 6*(Mcoarse(j,j) - Mcoarse(i,j));
    end
end


end

