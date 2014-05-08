function t = ext_knotvect(k,a,b,tau,M,ind )
L=length(tau)+1;
if ind==0
    %apro vettore esteso nodi
    t=a*ones(1,k);
    for i=1:L-1
        t=[t,tau(i)*ones(1,M(i))];
    end
    t=[t,b*ones(1,k)];
else
    %vettore esteso nodi chiuso
    t=a;
    for i=1:L-1
        t=[t,tau(i)*ones(1,M(i))];
    end
    t=[t,b];
    t=[zeros(1,k-1),t];
    idim=length(t);
    %ausiliari sinistri
    for i=k-1:-1:1
        t(i)=t(i+1)-(t(idim+i-k+1)-t(idim+i-k));
    end
    %ausiliari destri
    for i=1:k-1
        t(idim+i)=t(idim+i-1)+(t(i+k)-t(i+k-1));
    end
end
        



end

