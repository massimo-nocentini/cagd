%t vettore esteso nodi
%Q matrice punti di controllo
%ttab vettore di tabulazione
function C =deboor(t,Q,ttab )
[d,n]=size(Q);
t=t(:); t=t';
it=length(t);
k=it-n;
ttab=ttab(:); ttab=ttab';

ntab=length(ttab);
C=zeros(d,ntab);
ind=0;
for r=k:n
    if r<n
        isd=ttab>=t(r) & ttab<t(r+1);
    else
        isd=ttab>=t(r) & ttab<=t(r+1);
    end
    nloc=sum(isd);
    if nloc>0
        tloc=ttab(isd==1);
        Qloc=zeros(k,nloc,d);
        for i=1:d
            Qloc(:,:,i)=Q(i,r-k+1:r)'*ones(1,nloc);
        end
        for j=1:k-1
            alfa=zeros(k-1,nloc);
            for i=1:k-j
                alfa(i,:)=(tloc-t(i+r-k+j))/(t(i+r)-t(i+r-k+j));
                for s=1:d
                    Qloc(i,:,s)=(1-alfa(i,:)).*Qloc(i,:,s)+alfa(i,:).*Qloc(i+1,:,s);
                end
            end
        end
        for i=1:d
            C(i,ind+1:ind+nloc)=Qloc(1,1:nloc,i);
        end
        ind=ind+nloc;
    end
end



end

