function [ valid ] = check_schedule( Malap )

%Malap contains the delayed time steps of tiled algorithms

[rows, cols] = size(Malap);

ndiag = min(cols,rows);
%contains piv(i)
Eliminations = zeros(1,ndiag);

for j=1:ndiag,
    Eliminations(j) = j+1;
end

valid = 1;

for j=ndiag:-1:1,
    opcount = 0;
    
    if(rows==ndiag && j==ndiag)
        continue;
    end
    
    i = j;
    
    while(Eliminations(j) ~= rows+1),
        
        i=Eliminations(j);
        %find the number of operations at this timestep
        curval = Malap(i,j);
        k=i+1;
        
        while k<=rows,
            if(Malap(k,j)~=curval),
                break;
            end
            
            k=k+1;
        end
        
        opcount = k-i;
        %update piv(i)
        Eliminations(j) = k;
        
        
        %remove the elim
        for kk=k-1:-1:k-opcount,
            %            fprintf('Eliminating (%d,%d) with (%d,%d)\n',kk,j,kk-opcount,j);
            
            if(Malap(kk-opcount,j)>=Malap(kk,j)),
                Malap(kk-opcount,j) = Malap(kk,j);
                
                
                %count the updates of the elim
                for jj=j+1:cols-1,
                    if(kk>=jj),
                        if(Malap(kk,jj)-6>=Malap(kk,j)),
                            Malap(kk,jj) = Malap(kk,jj) -6;
                        else
                            valid = 0;
                            return;
                        end
                    end
                    
                    if(kk-opcount>=jj),
                        if(Malap(kk-opcount,jj)-6>=Malap(kk-opcount,j)),
                            Malap(kk-opcount,jj) = Malap(kk-opcount,jj) -6;
                        else
                            valid = 0;
                            return;
                        end
                    end
                    
                end
                %remove the elim
                Malap(kk,j) = Malap(kk,j)- 2;
                Malap(kk-opcount,j) = Malap(kk-opcount,j)- 2;
                
                
            else
                valid = 0;
                return;
            end
            
            
        end
        
    end
    
    for kk=j:rows,
        %count the updates
        for jj=j+1:cols-1,
            if(kk>=jj),
                if(Malap(kk,jj)-6>=Malap(kk,j)),
                    Malap(kk,jj) = Malap(kk,jj) -6;
                else
                    valid = 0;
                    return;
                end
            end
        end
        
        %remove the facto of the previous eliminations
        Malap(kk,j) = Malap(kk,j)- 4;
    end
end

end

