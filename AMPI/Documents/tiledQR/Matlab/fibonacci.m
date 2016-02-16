function [ z ] = fibonacci( p,q )
%MODI Summary of this function goes here
%   Detailed explanation goes here

%compute x least integer such that x(x+1)/2 >=  p-1
x=1;
while (x*(x+1)/2 < p-1),
    x=x+1;
end

z = zeros(p,q);

%fill first column
row=2;
for step = 1:x,
   for r =row:min(row + step -1,p),
        z(r,1) = x - step +1;
   end
   row = row + step;
end
%z(1,1) = z(2,1);

%complete other columns
for j=2:q,
    for i=j+1:p,
        z(i,j) = 2 + z(i-1,j-1);
    end
end

end

