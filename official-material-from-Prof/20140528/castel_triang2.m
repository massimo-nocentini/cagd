%
% function che esegue l'algoritmo di de Casteljau per valutare
% (contemporaneamente) in un vettore di tabulazione di coordinate
% baricentriche un patch di Bezier triangolare
%
% Input:
% -----
%
% n = grado del patch
% b = matrice dx[(n+1)(n+2)/2] contenente tutti i punti di controllo
%     di Bezier con punti dominio corrispondenti ordinati parallelamente 
%     al lato AB, da sinistra a destra, da AB verso C;
% ntab= numero di sottointervalli per la valutazione su ogni lato 
%       del dominio parametrico;
% V   = matrice 2x3 contenente le coordinate cartesiane (x e y) dei
%       vertici del dominio (solo nel caso funzionale)
% Output:
% ------
%
%  X = matrice 3x[(ntab+1)(ntab+2)/2] contenente le coordinate cartesiane
%      dei punti valutati sul Bezier patch 
%
%  tri = matrice ntriangx3 (ntriang=ntab^2) contenente la struttura dei triangoli per il patch 
%  U   = matrice 3 x ntab contnente le coordinate baricentriche dei punti del dominio parametrico di
%        valutazione (serve solo per il caso funzionale in cui i punti di
%        controllo diventano coefficienti di controllo
%
% trib = matrice ncx 3 (nc=n^2) contenente la struttura dei triangoli per
%        il net di Bezier (solo nel caso funzionale)
%
%  b   = matrice 3xntab contenente le coordinate cartesiane dei punti di
%        controllo (solo nel caso funzionale)
%
%   function [X,tri,U,trib,b]=castel_triang2(n,b,ntab,V)
%
  function [X,tri,U,trib,b]=castel_triang2(n,b,ntab,V)
  ir=1;
  [trib,Ub]=Ubar(n);
  
  if nargin==4 %caso funzionale
      ir=3;
      b(3,:)=b;
      b(1,:) = Ub(1,:).*V(1,1)+Ub(2,:).*V(1,2)+Ub(3,:)*V(1,3);
      b(2,:) = Ub(1,:).*V(2,1)+Ub(2,:).*V(2,2)+Ub(3,:)*V(2,3);
  end
  
  [d,ntot]=size(b);
  dims=(n+1)*(n+2)/2;
  
  if ntot> dims || ntot<dims
      error('numero di punti di controllo non adeguato a n')
  end
%
%  costruzione della matrice U (di dim. 3x[(ntab+1)(ntab+2)/2 ]) 
%  contenente le coordinate baricentriche dei punti di valutazione
%  e della corrispondente matrice di struttura tri
  [tri,U]=Ubar(ntab);
   N=size(U,2);
   X=zeros(3,N,ntot);
%
%  inizializzazione (in ogni punto di valutazione si inizializza 
%  con i punti di controllo)
%   
for id=ir:d
    X(id,:,:)=ones(N,1)*b(id,1:ntot);  % matrice dxNxntot
end 
for r=1:n  %passi di de Casteljau
    Xc=X; 
    nr=n-r+2;
    %    
    for k=0:n-r
        nrk=nr-k;
        for i=0:n-r-k  %percorrendoli con il solito ordinamento si calcolano
            % tutti i nuovi punti al passo r 
            ind1=sum(nrk+1:nr)+i+1;
            ind2=ind1+1;
            ind3=sum(nrk:nr)+i+1; % ind1,ind2,ind3 danno le posizioni
            % dei vecchi punti da usare         
            
            ind=sum(nrk:nr-1)+i+1;  % ind da' la posizione del nuovo punto
            for id=ir:d 
                X(id,:,ind)= U(1,:).*Xc(id,:,ind1)+U(2,:).*Xc(id,:,ind2)+...
                             U(3,:).*Xc(id,:,ind3);
            end
        end
    end
    %
end 
Xc=X; clear X;    

X(:,:)=squeeze(Xc(:,:,1)); % mi interessa solo l'ultimo valore
if nargin==4 %caso funzionale
    X(1,:)= U(1,:)*V(1,1)+U(2,:)*V(1,2)+U(3,:)*V(1,3);
    X(2,:)= U(1,:)*V(2,1)+U(2,:)*V(2,2)+U(3,:)*V(2,3);
end
%   
%

  