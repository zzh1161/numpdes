t = [0:0.01:4];
fi2 = tiledlayout(4,4,'TileSpacing','Compact','Padding','Compact');

f11 = (t-1).^2.*(t<=1);
f21 = (t-2).^2.*(t<=2);
f31 = (t-3).^2.*(t<=3);
f41 = (t-4).^2.*(t<=4);
f22 = (-2*t+3).*(t<=1)+(t-2).^2.*(t>1&t<=2);
f32 = (-2*t+5).*(t<=2)+(t-3).^2.*(t>2&t<=3);
f42 = (-2*t+7).*(t<=3)+(t-4).^2.*(t>3&t<=4);
f33 = (2).*(t<=1)+(-t.^2+2*t+1).*(t>1&t<=2)+(t-3).^2.*(t>2&t<=3);
f43 = (2).*(t<=2)+(-t.^2+4*t-2).*(t>2&t<=3)+(t-4).^2.*(t>3&t<=4);
f44 = (t-1).^2.*(t>1&t<=2)+(-2*t.^2+10*t-11).*(t>2&t<=3)+(t-4).^2.*(t>3&t<=4);

nexttile(1)
plot(t,f11);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

nexttile(5)
plot(t,f21);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

nexttile(6)
plot(t,f22);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

nexttile(9)
plot(t,f31);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

nexttile(10)
plot(t,f32);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

nexttile(11)
plot(t,f33);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

nexttile(13)
plot(t,f41);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

nexttile(14)
plot(t,f42);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

nexttile(15)
plot(t,f43);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

nexttile(16)
plot(t,f44);
axis equal;
axis([0 4 0 2.2]);
set(gca,'xtick',[0:1:4]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xticklabel',{'','$t_{i-1}$','$t_i$','$t_{i+1}$','$t_{i+2}$'});
set(gca,'yticklabel','');

title(fi2,'For n=2');