close all;
k=input('assegna ordine');
a=input('estremo sinistro');
b=input('estremo destro');
tau=input('vet nodi interni');
M=input('vet molteplicita');
ind=input('curva aperta(0) o chiusa');
t=ext_knotvect(k,a,b,tau,M,ind);
n=k+sum(M);
if ind==0
    disp(['assegna', num2str(n),'punti controllo (premi invio)']);
    pause
    axis;
    Q=ginput(n);
    close
    Q=Q';
else
    disp(['assegna', num2str(n-k+1),'punti controllo (premi invio)']);
     pause
    axis;
    Q=ginput(n-k+1);
    close
    Q=Q';
    Q=[Q,Q(:,1:k-1)];
end
ntab=input('tab');
ttab=linspace(a,b,ntab);
C=deboor(t,Q,ttab);
plot(C(1,:), C(2,:), 'b', Q(1,:), Q(2,:), 'r', Q(1,:), Q(2,:), 'ro');