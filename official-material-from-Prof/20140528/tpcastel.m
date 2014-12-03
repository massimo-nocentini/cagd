%
% Algoritmo di de Casteljau per le superfici tensor-product
% Si utilizza l'interpolazione bilineare fino a che e' possibile
% X in output e' di dimensioni mtabxntabx3 ma temporaneamente
% ha 5 dimensioni
%
% Input:
% -----
% b = matrice mxnx3 contenente le coordinate del net dei punti di controllo
% utab= vettore di tabulazione di mtab componenti in [0,1] del parametro u
% vtab= vettore di tabulazione di ntab componenti in [0,1] del parametro v
% ifunc = indicatore caso funzionale(1) o parametrico (0) 
%
% Nel caso funzionale ascisse e ordinate dei punti di controllo devono
% essere equispaziate
%
% Output:
% -------
% X = matrice mtabxntabx3 contenente le coordinate dei punti tabulati sulla
% superficie
%
% function X=tpcastel(b,utab,vtab,ifunc);
%
  function X=tpcastel(b,utab,vtab,ifunc)
%
  [m,n,d]=size(b);
  k=min(m,n);
  mtab=length(utab);
  ntab=length(vtab);
  [Vtab,Utab]=meshgrid(vtab,utab);  % dimensione mtabxntab
  Utm=1-Utab;
  Vtm=1-Vtab;
 
   X=zeros(mtab,ntab,d,m,n);
   
   ind=1;
   if ifunc==1
       ind=3;  % caso funzionale
   end
%
% inizializzazione
%
  for it=1:mtab
    for jt=1:ntab
      for id=ind:d
        X(it,jt,id,1:m,1:n)= squeeze(b(:,:,id));
      end 
    end
  end     
   
%
% si eseguono k-1 passi dell'algoritmo di de Casteljau
%  
 
  for r=1:k-1
     for i=1:m-r
       for j=1:n-r 
        for id=ind:d
          X(:,:,id,i,j)=     Utm .*Vtm .*squeeze(X(:,:,id,i,j))   +...
                           Utm .*Vtab.*squeeze(X(:,:,id,i,j+1)) + ...
                           Utab.*Vtm .*squeeze(X(:,:,id,i+1,j)) + ...
                           Utab.*Vtab.*squeeze(X(:,:,id,i+1,j+1));
         end                                         
       end
     end   
  end
%
%
  if m>n  %  occorrono altri passi in direzione u
     j=1;
     for r=1:m-k 
       for i=1:m-k+1-r
        for id=ind:d
          X(:,:,id,i,j)= Utm .*squeeze(X(:,:,id,i  ,j)) +...
                        Utab.*squeeze(X(:,:,id,i+1,j)); 
        end  
       end
     end 

   elseif n>m   %  occorrono altri passi in direzione v  
    i=1;
     for r=1:n-k 
       for j=1:n-k+1-r
         for id=ind:d
          X(:,:,id,i,j)= Vtm .*squeeze(X(:,:,id,i,j  )) +...
                        Vtab.*squeeze(X(:,:,id,i,j+1));  
         end             
       end
     end  
   else  % se m=n ho concluso
  end 
       
  if ifunc==1 % ascisse e ordinate dei punti di controllo equispaziate
      X(:,:,1,1,1)= (1-Utab)*b(1,1,1)+Utab*b(1,end,1);
      X(:,:,2,1,1)= (1-Vtab)*b(1,1,2)+Vtab*b(end,1,2);
  end
   
  X=squeeze(X(:,:,:,1,1));     
  
 
    
  end
 