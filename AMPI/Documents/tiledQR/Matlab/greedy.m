function [ M ] = greedy( p,q )
   nZ = zeros(1,q+1);
   nZ(1) = p;
   M = zeros(p,q);

   if (p==q), q = q-1; end

   time = 1;

   while ( nZ(q+1) < p-q );
   
      for j = q:-1:1,
         z = floor( ( nZ(j) - nZ(j+1) ) / 2 ) ;
         M(p-nZ(j+1)-z+1:p-nZ(j+1), j) = time;
         nZ(j+1) = nZ(j+1) + z ;
      end

      time = time+1;

   end

end
