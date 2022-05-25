%corollary 4.58
x = [-8:0.02:8];

y1 = ( 2.29911047382 + 0.862166427682*x.^1 + 0.10777080346*x.^2 + 0.00449045014418*x.^3 ).*(x>=-8 & x<-7) ...
+( -1.7566527754 + -0.876017821981*x.^1 + -0.140541232206*x.^2 + -0.00733393250659*x.^3 ).*(x>=-7 & x<-6) ...
+( 0.795267511986 + 0.399942321709*x.^1 + 0.0721187917425*x.^2 + 0.00448051326832*x.^3 ).*(x>=-6 & x<-5) ...
+( 0.317453206474 + 0.113253738403*x.^1 + 0.0147810750811*x.^2 + 0.000657998824229*x.^3 ).*(x>=-5 & x<-4) ...
+( 0.768288486936 + 0.451380198749*x.^1 + 0.0993126901677*x.^2 + 0.00770230008144*x.^3 ).*(x>=-4 & x<-3) ...
+( 0.736956351579 + 0.420048063392*x.^1 + 0.0888686450487*x.^2 + 0.00654185062377*x.^3 ).*(x>=-3 & x<-2) ...
+( 1.54307569068 + 1.62922707205*x.^1 + 0.693458149376*x.^2 + 0.107306768012*x.^3 ).*(x>=-2 & x<-1) ...
+( 1 + -5.55111512313e-17*x.^1 + -0.935768922671*x.^2 + -0.435768922671*x.^3 ).*(x>=-1 & x<0) ...
+( 1 + -5.55111512313e-17*x.^1 + -0.935768922671*x.^2 + 0.435768922671*x.^3 ).*(x>=0 & x<1) ...
+( 1.54307569068 + -1.62922707205*x.^1 + 0.693458149376*x.^2 + -0.107306768012*x.^3 ).*(x>=1 & x<2) ...
+( 0.736956351579 + -0.420048063392*x.^1 + 0.0888686450487*x.^2 + -0.00654185062377*x.^3 ).*(x>=2 & x<3) ...
+( 0.768288486936 + -0.451380198749*x.^1 + 0.0993126901677*x.^2 + -0.00770230008144*x.^3 ).*(x>=3 & x<4) ...
+( 0.317453206474 + -0.113253738403*x.^1 + 0.0147810750811*x.^2 + -0.000657998824229*x.^3 ).*(x>=4 & x<5) ...
+( 0.795267511986 + -0.399942321709*x.^1 + 0.0721187917425*x.^2 + -0.00448051326832*x.^3 ).*(x>=5 & x<6) ...
+( -1.7566527754 + 0.876017821981*x.^1 + -0.140541232206*x.^2 + 0.00733393250659*x.^3 ).*(x>=6 & x<7) ...
+( 2.29911047382 + -0.862166427682*x.^1 + 0.10777080346*x.^2 + -0.00449045014418*x.^3 ).*(x>=7 & x<8);


%corollary 4.59
y2 = ( 0.764964632258 + 0.218561323502*x.^1 + 0.0156115231073*x.^2 ).*(x>=-7 & x<-6) ...
+( -0.0984793105793 + -0.0692533241102*x.^1 + -0.0083730308604*x.^2 ).*(x>=-6 & x<-5) ...
+( 0.246725743027 + 0.0688286973323*x.^1 + 0.00543517128384*x.^2 ).*(x>=-5 & x<-4) ...
+( 0.442920784453 + 0.166926218045*x.^1 + 0.017697361373*x.^2 ).*(x>=-4 & x<-3) ...
+( 0.504743100688 + 0.208141095535*x.^1 + 0.0245665076213*x.^2 ).*(x>=-3 & x<-2) ...
+( 1.46292243124 + 1.16632042609*x.^1 + 0.264111340259*x.^2 ).*(x>=-2 & x<-1) ...
+( 0.879762218196 + 0*x.^1 + -0.319048872784*x.^2 ).*(x>=-1 & x<0) ...
+( 0.879762218196 + 0*x.^1 + -0.319048872784*x.^2 ).*(x>=0 & x<1) ...
+( 1.46292243124 + -1.16632042609*x.^1 + 0.264111340259*x.^2 ).*(x>=1 & x<2) ...
+( 0.504743100688 + -0.208141095535*x.^1 + 0.0245665076213*x.^2 ).*(x>=2 & x<3) ...
+( 0.442920784453 + -0.166926218045*x.^1 + 0.017697361373*x.^2 ).*(x>=3 & x<4) ...
+( 0.246725743027 + -0.0688286973323*x.^1 + 0.00543517128384*x.^2 ).*(x>=4 & x<5) ...
+( -0.0984793105793 + 0.0692533241102*x.^1 + -0.0083730308604*x.^2 ).*(x>=5 & x<6) ...
+( 0.764964632258 + -0.218561323502*x.^1 + 0.0156115231073*x.^2 ).*(x>=6 & x<7);

origin = 1./(1+x.^2);

plot(x,origin,x,y1,x,y2);
axis([-10 10 0 1]);
legend('origin','Corollary4.58','Corollary4.59');
title("ASSIGNMENT.C");

