warning off
clear 
clc
global Lt La1 La T;
global f fx fu L Lx Lu cL cLx cLu n lx lu Nh NO;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NO=1;%NO+1为有限单元上取的配置点数,配置点适当增加，则求解越精确，但是求解时间越长
Nh=1;%有限元长度，有限单元数过大会导致微分方程模拟失败

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=0.01;beta=4e-5;rho=0.00;%参数设置
% f=@(x,u)[-beta*x(1)*x(2)-u*x(1),beta*x(1)*x(2)-mu*x(2),-rho*x(3)];
% fx=@(x,u)[-beta*x(2)-u,-beta*x(1),0;beta*x(2),beta*x(1)-mu,0;0,0,-rho];
% fu=@(x,u)[-x(1);0;0];
f=@(x,u)[-beta*(1-u(1))*x(1)*x(2)-u(2)*x(1),beta*(1-u(1))*x(1)*x(2)-mu*x(2),-rho*x(3)];
fx=@(x,u)[-beta*(1-u(1))*x(2)-u(2),-beta*(1-u(1))*x(1),0;beta*(1-u(1))*x(2),beta*(1-u(1))*x(1)-mu,0;0,0,-rho];
fu=@(x,u)[beta*x(1)*x(2) -x(1);-beta*x(1)*x(2) 0;0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=1;ku=1;%被积函数参数设置
%  L=@(x,u)(c*x(2)+ku*u)*x(3);Lx=@(x,u)[0,c*x(3),(c*x(2)+ku*u)];Lu=@(x,u)ku*x(3);%目标中的积分函数
 L=@(x,u)c*x(2)*x(3)+ku*u(2)*x(3);Lx=@(x,u)[0,c*x(3),(c*x(2)+ku*u(2))];Lu=@(x,u)[0,ku*x(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cL cLx cLu;
cmax=10;%单日最大隔离人数
% cL=@(x,u)u*x(1)-cmax;cLx=@(x,u)[u,0,0];cLu=@(x,u)x(1);%逐点约束函数及其对x,u的梯度函数
cL=@(x,u)u(2)*x(1)-cmax;cLx=@(x,u)[u(2),0,0];cLu=@(x,u)[0,x(1)];%逐点约束函数及其对x,u的梯度函数
cL=@(x,u)[u(2)*x(1)-cmax;u(2)*x(1)-cmax];cLx=@(x,u)[u(2),0,0;u(2),0,0];cLu=@(x,u)[0,x(1);0,x(1)];%逐点约束函数及其对x,u的梯度函数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0=[1000,100,1];
n=10;%有限单元个数，每一个有限单元上取相同的控制，因此n也相当于控制分段数，
lu=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=0.0001*rand(n,lu);%随机产生一个可行的初始控制。
Ul=0*ones(size(U));%控制下界
UM=0.2*ones(size(U));%控制上界为0.2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%求解算法设置，
%可选的非线性约束优化算法有内点法'interior-point',序贯二次规划法'sqp',信赖域法'trust-region-reflective',有效集法'active-set'
algorithm='active-set';algorithm='interior-point';
%优化迭代过程显示控制
% @optimplotx                    plots the current point
% @optimplotfunccount            plots the function count
% @optimplotfval                 plots the function value
% @optimplotconstrviolation      plots the maximum constraint violation
% @optimplotstepsize             plots the step size
% @optimplotfirstorderopt        plots the first-order optimality measure
plotfcns={@optimplotfval,@optimplotconstrviolation,@optimplotfirstorderopt,@optimplotx};
%最大调用次数
maxfunevals=3000;

disp('程序输入完成，转入input_check.m进行输入检查')
pause(0.5)
input_check
disp('检查通过,转入ultimate.m进行最优控制的求解')
pause(0.5)
ultimate
    


