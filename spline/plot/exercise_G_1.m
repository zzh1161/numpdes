x = [0:0.02:3];
f1_1 = (-x+1).*(x<=1);
f2_1 = (-x+2).*(x<=2);
f3_1 = (-x+3).*(x<=3);
f1_2 = (-x+2).*(x>=1 & x<2)+(1).*(x<1);
f2_2 = (-x+3).*(x>=2 & x<3)+(1).*(x<2);
f1_3 = (x-1).*(x>=1 & x<2)+(-x+3).*(x>=2 & x<3);

fi = tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');

nexttile
plot(x,f1_1);
axis equal;
axis([0 3 0 1]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$'});
set(gca,'yticklabel','');

nexttile(4)
plot(x,f2_1);
axis equal;
axis([0 3 0 1]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$'});
set(gca,'yticklabel','');

nexttile(5)
plot(x,f1_2);
axis equal;
axis([0 3 0 1]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$'});
set(gca,'yticklabel','');

nexttile(7)
plot(x,f3_1);
axis equal;
axis([0 3 0 1]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$'});
set(gca,'yticklabel','');

nexttile(8)
plot(x,f2_2);
axis equal;
axis([0 3 0 1]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$'});
set(gca,'yticklabel','');

nexttile(9)
plot(x,f1_3);
axis equal;
axis([0 3 0 1]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$'});
set(gca,'yticklabel','');

title(fi,'For n=1');