Mesh-networks

To measure the Gateway- average path length with NLL (No. of long ranged links) varied for non- uniform grid network

Multi-Gateway Aware LL addition Strategy (M-GAS)

-Intra-region LLs:

For any nodes i and j lying in same region, check following conditions:

deucl(i,j) lies between R1 and R2. (R1=short-link transmission range; R2=maximum length of long-ranged links)

Each node (i and j) has an unoccupied LL radio.

|d(i)-d(j)|>=delta(h) ; d(i) = shortest path length between node i and their closest Gateway (in hops);h is input parameter.

-Inter-region LLs:

For any nodes i and j lying in different regions, check following conditions:

deucl(i,j) lies between R1 and R2. (R1=short-link transmission range; R2=input parameter)

Each node (i and j) has an unoccupied LL radio.

|d(i)-d(j)|>=delta(h) d(i) = shortest path length between node i and its closest Gateway G(i) (wrt. number of hops) in hops; ?h is input parameter.
