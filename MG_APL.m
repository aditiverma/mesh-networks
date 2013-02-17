
function [aver_G_APL] = MG_APL(A3)

 N=size(A3,2);
 D=A3;

 D(find(D==0))=inf;  
 for i=1:N
     D(i,i)=0;       
 end   
 
 
 %calculating shortest path between nodes and gateways 
 for k=1:N
     for i=1:N
         for j=1:N
             if D(i,j)>D(i,k)+D(k,j)
                D(i,j)=D(i,k)+D(k,j);
             end
         end
     end
 end
 
 for j=1:N-4
     node_gateway_min(j) = min(D(j,N-3:N));
     aver_G_APL = mean(node_gateway_min);
 
 end


