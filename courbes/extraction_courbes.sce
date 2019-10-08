x=fscanfMat('z.txt') ; y=fscanfMat('lambda.txt') ; P=fscanfMat('Psdbm.txt') ;
x2 = [0:0.2:17.6] ; y2=[1.52:0.001:1.57] ;
P2 = P(:,281:4:481) ;

fenetre = scf(100001);
clf(fenetre,"reset");

Titre  = "Directions d''arrivée des signaux (k.d.sin(teta))";
my_handle.figure_name = Titre;

xlabel('Lambda (nm)') ; ylabel('Ps (dBm)');
plot(y2,P2(5,:),y2,P2(10,:),y2,P2(15,:),y2,P2(20,:),y2,P2(25,:));
xgrid();
legend('1m','2m','3m','4m','5m');
