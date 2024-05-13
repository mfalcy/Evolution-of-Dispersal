library(igraph)
edges<-c(1,2,2,3,2,4,3,5,3,6,4,7,4,8,5,9,6,10,7,11,8,12,9,13,9,14,10,15,10,16,11,17,11,18,12,19,12,20)#All upstream pairwise connections (1 to 2, 2 to 3, 2 to 4, etc.)
g<-graph(edges)
plot(g,edge.arrow.size=0.3)
#basic calls
all_simple_paths(g,from=1,to=15)#Returns all the nodes between node 1 and node 15
unlist(adjacent_vertices(g,3,mode="out"))#upstream choices from a given node.

#In the simulation below, an animal is from node 18 (home=18, Line 14). We watch where it goes over 10000 replicates.  
Nreps<-10000
out<-vector()#will hold simulation results of Nreps of the returns to home node 18
p<-0.85#probability of correctly choosing Left vs. Right
for (rep in 1:Nreps){
home<-18
path<-unlist(all_simple_paths(g,from=1,to=home))
spawn<-c(13,14,15,16,17,18,19,20)#all the spawning nodes
pathi<-2 #start at node 2.
while(length(intersect(pathi,spawn))<1){
    #on correct path?
    if (sum(pathi==path)>0){
     #is there a fork?
      if (length(unlist(adjacent_vertices(g,pathi,mode="out")))>1){
        true<-intersect(unlist(adjacent_vertices(g,pathi,mode="out")),path)
        options<-unlist(adjacent_vertices(g,pathi,mode="out"))
        pathi<-options[true!=options]#default to wrong choice
        if (rbinom(1,1,p)==1){
          pathi<-options[true==options]#fix choice
        }
      }
      #no fork, so no real choice
      else {
        pathi<-unlist(adjacent_vertices(g,pathi,mode="out"))
      }
    }
    #not on the right path
    else{
      #is there a fork?
      if (length(unlist(adjacent_vertices(g,pathi,mode="out")))>1){
        pathi<- unlist(adjacent_vertices(g,pathi,mode="out"))[rbinom(1,1,0.5)+1]#chooses randomly among two upstream nodes
      }
      #no fork, so no real choice 
      else{
        pathi<-unlist(adjacent_vertices(g,pathi,mode="out"))
      }
    }
}
out[rep]<-pathi
}
hist(out,breaks=seq(12.5,20.5))
  
  