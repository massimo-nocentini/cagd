%
% function che calcola un set uniforme di punti in un triangolo di
% riferimento espressi in coordinate baricentriche 
%
% Input:
% -----
% 
% ntab= numero di sottointervalli per la valutazione su ogni lato 
%       del dominio parametrico;
% 
% Output:
% ------
%
%  % tri = matrice con tante righe quanti sono i triangoli (ntab^2) e 3 colonne 
%      (specifica la struttura della triangolazione uniforme) 
%
%  U  = matrice 3x[(ntab+1)(ntab+2)/2] contenente le coordinate
%       baricentriche dei punti uniformememte distribuiti sul triangolo 
%
%   function [tri,U]=Ubar(ntab)
%
  function [tri,U]=Ubar(ntab)
   
tri = zeros(ntab^2,3);

count=0;
for kt=0:ntab  
    np=ntab-kt+1;
    U(1:3,count+1:count+np)= ...
        [np-1:-1:0;0:np-1;kt*ones(1,np)]/ntab;    
    count=count+np;
end  
N=size(U,2);
%
%   matrice di struttura della triangolazione (per la trisurf)
%
count=0;
for kt=0:ntab-2  
    nk = ntab+2-kt;
    sm = sum(nk:ntab+1);
    for it=0:ntab-kt-2
        ind  = sm+it+1;
        count= count+1;       
        tri(count,:)=[ind,ind+1,ind+nk-1];
        count=count+1;       
        tri(count,:)=[ind+1,ind+nk-1,ind+nk];
    end
    count=count+1;       
    tri(count,:)=[ind+1,ind+2,ind+nk];
end
count=count+1;    
tri(count,:)=[N-2,N-1,N];     
   