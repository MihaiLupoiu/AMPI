function [ minM ] = quadridiagonal( p,q )

%generate every valid permutations

validperm= [];
for i = 2:1:4,
    for j = 2:1:4,
        for k = 2:1:4,
            if(j~=i && k~=j && i~=k)
                if(k<j)
                    scheme = [k 1 j 1 i 1];
                    validperm = [validperm ; scheme];
                end
                
                if(k>i)
                    scheme = [k i j 1 i 1];
                    validperm = [validperm ; scheme];
                    if(k>j && j>i)
                        scheme = [k i j i i 1];
                        validperm = [validperm ; scheme];
                    end
                end
                if(k>j)
                    scheme = [k j j 1 i 1];
                    validperm = [validperm ; scheme];
                end
                
            end
        end
    end
end


%
eliminations = pick(1:10,min(q,p),'r');
[rows,cols]=size(eliminations);


mindiag = inf;
minexpe = -1;

minM = zeros(p,q);


for expe = 1:1:rows,
    
    permutM = zeros(p,q);
    q = zeros(p,q);
    
    for i = 1:1:min(q,p),
        %for j = i:1:p,
        for k = i:1:i+3,
            if(k<=min(q,p))
                %factor
                q(k,i) = q(k,i)+ 4;
                %update
                for jj = i+1:1:q,
                    q(k,jj) = max(q(k,jj),q(k,i)) + 6;
                end
            end
        end
        
        %eliminate
        schemeindex = eliminations(expe,i);
        %schemeindex=1;
        scheme = validperm(schemeindex,:);
        for op = 1:2:6,
            
            bl = i - 1 + scheme(op);
            tl = i - 1 + scheme(op+1);
            if(bl<=p && tl<=p)
                
                maxval = max(q(tl,i),q(bl,i));
                q(tl,i) = maxval + 2;
                q(bl,i) = maxval + 2;
                
          
                %update
                for jj = i+1:1:q,
                    maxval = max(q(tl,jj),q(bl,jj));
                    q(tl,jj) = max(max(q(tl,i),q(bl,i)),maxval) + 6;
                    q(bl,jj) = max(max(q(tl,i),q(bl,i)),maxval) + 6;
                end
            end
        end
        
        %end
        %break
    end
    
        for i = 1:1:p,
                for jj = i+1:1:q,
                   q(i,jj) = 0;
                end
        end
        
    %if (diagval<=mindiag)
    %    mindiag = diagval;
    %    minexpe = expe;
    %    minM = q;
    %end
    

    
    if (q(min(q,p),min(q,p))<=mindiag)
        mindiag = q(min(q,p),min(q,p));
        minexpe = expe;
        minM = q;
%        for i = 1:1:min(q,p),
%            schemeindex = eliminations(minexpe,i);
%            scheme = validperm(schemeindex,:)
%        end
%        mindiag
%        minM
    end

    
end


for i = 1:1:min(q,p),
    schemeindex = eliminations(minexpe,i);
    scheme = validperm(schemeindex,:)
    %schemeindex=1;
    %scheme = validperm(schemeindex,:);
end
mindiag
minM
for i = min(q,p):-1:2,
    minM(i,i) = minM(i,i) - minM(i-1,i-1);
end
minM
end