

%--------------------
% Code details: To vary NLL in a WMN with four gateways with proper load balancing.
% Inputs will be asked from the user.
% --------------------
clear;
clear all;
N=input('Enter N (Note that number of nodes = N*N): '); 
%N=10;
R1=input('Enter R1 (Range for short-range links): ');
%R1=1;
R2=input('Enter R2 (Max range for long-ranged links): ');
%R2=8;
max_num_LL_radios=input('Enter the max. no. of LL radios that a node can have:');
%max_num_LL_radios=2;
dmin=input('enter minimum difference between number of hops of two nodes from a gateway.');
%dmin=0;
delta_n=2;
%Assigning X and Y co-ordinates to each node FOR GRID
for j=1:N
    for i=1:N
        x(i+N*(j-1))=i; %assigning X-coordinate
        y(i+N*(j-1))=j; %assigning Y-coordinate
    end      %hence x[] and y[] are populated
end % GRID mode asignment of A ends      

%-------Defining the Adjacency matrix ----------
A=zeros(N*N); % N^2 by N^2 adjacency matrix, initialized to zeros
for i=1:(N*N)
    for j=i+1:(N*N) %'i'th node to all nodes numbered above it
        dist_sq = (x(i)-x(j))^2 + (y(i)-y(j))^2 ;
        if( ( dist_sq <= R1*R1) && (i~=j))
            A(i,j)=1;
            A(j,i)=1;
        end
    end
end
% 'A' has been filled appropriately, without long-range links

 a=(1+N)/2;
 b=(1+N)/2;

x_gateway(1)=(1+a)/2 ;
x_gateway(2)=(((N+a)/2 + N)/2+ N)/2;
x_gateway(3)=((1+a)/2 +a)/2 ;
x_gateway(4)=((N+a)/2+N)/2;
y_gateway(1)=(1+b)/2 ;
y_gateway(2)=((1+b)/2 +1)/2 ;
y_gateway(3)=((b+N)/2 +b)/2  ;
y_gateway(4)=(((b+N)/2 + N)/2 +N)/2 ; 



%adding gateway nodes in adjacency matrix A2
A2=A;
for i=1:N*N
    for k= 1:4
       d= (x(i) - x_gateway(k))^2 + (y(i) - y_gateway(k))^2; 
       if(d<=R1*R1)
       A2(N^2 + k,i)=1;
       A2(i,N^2 + k)=1;
       else
       A2(N^2 + k,i)=0;
       A2(i,N^2 + k)=0;
       end                   
    end
end       
           
%defining a matrix for no. of hops between nodes and gateways
%without long links
N1=size(A2,2);
D=A2;

D(find(D==0))=inf;  
for i=1:N1
    D(i,i)=0;       
end   

 %calculating shortest path between nodes and gateways without
 %LLs
for k=1:N1
   for i=1:N1
       for j=1:N1
           if D(i,j)>D(i,k)+D(k,j)
                 D(i,j)=D(i,k)+D(k,j);
           end
       end
   end
end
for p= 1:4
   g(:,p)= D(:,N*N +p); %g is used for checking condition of dmin in intra-node links
end
affliated_count= [0 0 0 0];  
 g2=g;
 for k2=1:N*N
   min_hops_for_k2th_node = min(g2(k2,:));
   location_min_hops = find(g2(k2,:)==min_hops_for_k2th_node);
   num_min_hops_for_k2th_node = size(location_min_hops,2) ;
   if (num_min_hops_for_k2th_node ==1) % the case when the minimum is unique
      index= find( g2(k2,:)==min_hops_for_k2th_node ); 
      affliation(k2) = index;
      affliated_count(index)=affliated_count(index)+1;
   end
   if(num_min_hops_for_k2th_node == 2) 
       affliation(k2) = location_min_hops(1);
       r1= rand();
       if (r1<0.5)
           affliation(k2) = location_min_hops(2);
       end
       affliated_count(affliation(k2))=affliated_count(affliation(k2))+1;
   end
   if(num_min_hops_for_k2th_node == 3)
       affliation(k2) = location_min_hops(1);
       r1= rand();
       if (r1<0.33)
           affliation(k2) = location_min_hops(2);
       end
       if (r1>=0.33 && r1<0.66)
           affliation(k2) = location_min_hops(3);
       end
       affliated_count(affliation(k2))=affliated_count(affliation(k2))+1;
   end
   if(num_min_hops_for_k2th_node == 4)
       affliation(k2) = location_min_hops(1);
       r1= rand();
       if (r1<0.25)
           affliation(k2)= location_min_hops(2);
       end
       if (r1>=0.25 && r1<0.5)
           affliation(k2) = location_min_hops(3);
       end
       if (r1>=0.5 && r1<0.75)
           affliation(k2) = location_min_hops(4);
       end
       affliated_count(affliation(k2))=affliated_count(affliation(k2))+1;
   end
 end % for k2=1:N*N

A_safe=A;  
           


links_intra_2=zeros(20,13);
links_inter_2=zeros(20,13);  
n_long_links=zeros(20,13);  
 
argument = 0;

for N_LL=0:8:(N^2)   %main outer loop        
%  N_LL=16;
    argument = argument +1;

    for repitition = 1:20    
        A=A_safe;

        %----Modifying A to introduce long-links in system by link-addition 
        B=A; %copying A to B. B will be used to store the connections which have been considered to be LL
    %     n_long_links(repitition,argument)=0; %initializing n_long_links = number of LLs added

        % -- link addition -- % `
        for i=1:N*N
            for j=1:i
                B(i,j)=1; % to ensure that the diagnoal items and left triangle ones don't remain '0'
            end
        end
    n_unconsidered_LL = length(find(B==0)); % these many links are remaining to be considered to be LL
        LL_array=randperm(n_unconsidered_LL ); % this array contains random sequence for LL-consideration
        [n1, n2] = find(B==0); % to get which nodes are connected by unconsidered links
        % n1 and n2 would be 2 arrays with corresponding (respective) node no.s
        % so that we can identify the connecting link


        m = N_LL/2; 
           C1 = zeros(N*N,N*N);
            A1=A;
            n=abs(N_LL-m); 
            for num_each_link_2=1:n_unconsidered_LL
                link_added=false;
                link_num_2 = LL_array(num_each_link_2);  
                  % For intra region LLs            
		    if (affliation(n1(link_num_2)) == affliation(n2(link_num_2)) && links_intra_2(repitition,argument) < m)             % change made. 
                    dis_sq=( x(n1(link_num_2))-x(n2(link_num_2)) )^2 + (y(n1(link_num_2))-y(n2(link_num_2)) )^2;
                   [value1 index1]=min(g(n1(link_num_2),:));
                    [value2 index2]=min(g(n2(link_num_2),:));
                    d1 = abs(g(n1(link_num_2),index1)-g(n2(link_num_2),index2));%d1 is difference in no. of hops of two nodes from a gateway
                    if ((length(find(C1(n1(link_num_2),:)==2))<max_num_LL_radios) &&(length(find(C1(n2(link_num_2),:)==2)) <max_num_LL_radios)&& (dis_sq>R1*R1) && (dis_sq<R2*R2) && d1>=2 )
                        A1(n1(link_num_2),n2(link_num_2)) =1;
                        A1(n2(link_num_2),n1(link_num_2)) =1;
                        link_added=true;
                           C1(n1(link_num_2),n2(link_num_2)) =2;
                          C1(n2(link_num_2),n1(link_num_2)) =2;
                        %disp(sprintf('Intra link added \n'));
                        n_long_links(repitition,argument)=n_long_links(repitition,argument)+1;
                        links_intra_2(repitition,argument)=links_intra_2(repitition,argument) + 1;
                    end %end of wireless constraints for inter-node links in same quarter

                    % NOW FOR INTER - LINKS
                elseif (affliation(n1(link_num_2)) ~= affliation(n2(link_num_2)) && links_inter_2(repitition,argument) < n)
                    %B1(n1(link_num_2),n2(link_num_2)) = 1; % indicating we have considered this link for LL -- although this will not be used
                    dis_sq=( x(n1(link_num_2))-x(n2(link_num_2)) )^2 + (y(n1(link_num_2))-y(n2(link_num_2)) )^2;
                    [value1 index1]=min(g(n1(link_num_2),:));
                    [value2 index2]=min(g(n2(link_num_2),:));

                    d1 = abs(g(n1(link_num_2),index1)-g(n2(link_num_2),index2) );%d1 is difference in no. of hops of two nodes from a gateway

                    if ( (length(find(C1(n1(link_num_2),:)==2))<max_num_LL_radios) &&(length(find(C1(n2(link_num_2),:)==2)) <max_num_LL_radios)&& (dis_sq>R1*R1) && (dis_sq<R2*R2) &&  (d1>=dmin) )
                        A1(n1(link_num_2),n2(link_num_2)) =1;
                        A1(n2(link_num_2),n1(link_num_2)) =1;
                        link_added=true;
                         C1(n1(link_num_2),n2(link_num_2)) =2;
                        C1(n2(link_num_2),n1(link_num_2)) =2; 
                        %disp(sprintf('Inter link added \n'));

                        n_long_links(repitition,argument)=n_long_links(repitition,argument)+1;
                        links_inter_2(repitition,argument)=links_inter_2(repitition,argument) + 1;
                    end %end of wireless constraints for inter-node links in different quarters.  
                end %end of elseif

                if ((links_inter_2(repitition,argument) ==n) && (links_intra_2(repitition,argument) == m))
                    break;
                end


                if(link_added) % update hopcounts in g, update affliations and update affliated count
                        A2=A1;
                        for i=1:N*N
                            for k= 1:4
                                d= (x(i) - x_gateway(k))^2 + (y(i) - y_gateway(k))^2; 
                                if(d<=R1*R1)
                                    A2(N^2 + k,i)=1;
                                    A2(i,N^2 + k)=1;
                                else
                                    A2(N^2 + k,i)=0;
                                    A2(i,N^2 + k)=0;
                                end                   
                            end
                        end       
                        %defining a matrix for no. of hops between nodes and gateways
                        %without long links
                        N1=size(A2,2);
                        D=A2;                
                        D(find(D==0))=inf;  
                        for i=1:N1
                            D(i,i)=0;       
                        end                   
                        %calculating shortest path between nodes and gateways without
                        %LLs
                        for k=1:N1
                            for i=1:N1
                                for j=1:N1
                                    if D(i,j)>D(i,k)+D(k,j)
                                        D(i,j)=D(i,k)+D(k,j);
                                    end
                                end
                            end
                        end
                        for p= 1:4
                            g(:,p)= D(:,N*N +p); %g is used for checking condition of dmin in intra-node links
                        end   
                        affliated_count= [0 0 0 0];  
                         g2=g; 
                         for k2=1:N*N
                           min_hops_for_k2th_node = min(g2(k2,:));
                           location_min_hops = find(g2(k2,:)==min_hops_for_k2th_node);
                           num_min_hops_for_k2th_node = size(location_min_hops,2) ;
                           if (num_min_hops_for_k2th_node ==1) % the case when the minimum is unique
                              index= find( g2(k2,:)==min_hops_for_k2th_node ); 
                              affliation(k2) = index;
                              affliated_count(index)=affliated_count(index)+1;
                           end
                           if(num_min_hops_for_k2th_node == 2) 
                               affliation(k2) = location_min_hops(1);
                               r1= rand();
                               if (r1<0.5)
                                   affliation(k2) = location_min_hops(2);
                               end
                               affliated_count(affliation(k2))=affliated_count(affliation(k2))+1;
                           end
                           if(num_min_hops_for_k2th_node == 3)
                               affliation(k2) = location_min_hops(1);
                               r1= rand();
                               if (r1<0.33)
                                   affliation(k2) = location_min_hops(2);
                               end
                               if (r1>=0.33 && r1<0.66)
                                   affliation(k2) = location_min_hops(3);
                               end
                               affliated_count(affliation(k2))=affliated_count(affliation(k2))+1;
                           end
                           if(num_min_hops_for_k2th_node == 4)
                               affliation(k2) = location_min_hops(1);
                               r1= rand();
                               if (r1<0.25)
                                   affliation(k2)= location_min_hops(2);
                               end
                               if (r1>=0.25 && r1<0.5)
                                   affliation(k2) = location_min_hops(3);
                               end
                               if (r1>=0.5 && r1<0.75)
                                   affliation(k2) = location_min_hops(4);
                               end
                               affliated_count(affliation(k2))=affliated_count(affliation(k2))+1;
                           end
                         end % for k2=1:N*N              
                end %end of if link_added            
            end %end for one unconsidered LL.     


            %adding gateway nodes in adjacency matrix
            for i=1:N*N
                for k= 1:4
                    d= (x(i) - x_gateway(k))^2 + (y(i) - y_gateway(k))^2; 
                    if(d<=R1*R1)
                        A1(N^2 + k,i)=1;
                        A1(i,N^2 + k)=1;
                    else
                        A1(N^2 + k,i)=0;
                        A1(i,N^2 + k)=0;
                    end                   
                end
            end
            [average_G_APL(repitition,argument)] = MG_APL(A1); 

             [P] = MG_APL2(A1);
             g2 = P(:,N*N+1:N*N+4) ; % g2 stores last 4 columns 
			 count= [0 0 0 0];  
             for k2=1:N*N
               min_hops_for_k2th_node = min(g2(k2,:));
               location_min_hops = find(g2(k2,:)==min_hops_for_k2th_node);
               num_min_hops_for_k2th_node = size(location_min_hops,2) ;
               if (num_min_hops_for_k2th_node ==1) % the case when the minimum is unique
                  index= find( g2(k2,:)==min_hops_for_k2th_node ); 
                  count(index)=count(index) + 1;
               end
               if(num_min_hops_for_k2th_node == 2) 
                   index = location_min_hops(1);
                   r1= rand();
                   if (r1<0.5)
                       index = location_min_hops(2);
                   end
                  count(index)=count(index) + 1;
               end
               if(num_min_hops_for_k2th_node == 3)
                   quart(k2) = location_min_hops(1);
                   r1= rand();
                   if (r1<0.33)
                       quart(k2) = location_min_hops(2);
                   end
                   if (r1>=0.33 && r1<0.66)
                       quart(k2) = location_min_hops(3);
                   end
                   count(index)=count(index) + 1;	
               end
               if(num_min_hops_for_k2th_node == 4)
                   quart(k2) = location_min_hops(1);
                   r1= rand();
                   if (r1<0.25)
                       quart(k2) = location_min_hops(2);
                   end
                   if (r1>=0.25 && r1<0.5)
                       quart(k2) = location_min_hops(3);
                   end
                   if (r1>=0.5 && r1<0.75)
                       quart(k2) = location_min_hops(4);
                   end
                   count(index)=count(index) + 1;	
               end
             end % for k2=1:N*N

               %calculating standard deviation
			sdev(repitition,argument)= std(count);
          end %for the repititions
end % for N_LL
     
     sd=mean(sdev);
   
%  avg_n_long_links(argument)=mean(n_long_links);         
%  avg_aver_D( argument ) = mean (aver_G_APL); % mean over 20 iterations         
% end of main outer loop
 disp('simulation done_21aug!');
 save NLLvaried_MGAS_nonuniform_10x10grid_equalinterintra_dmin0_21aug



