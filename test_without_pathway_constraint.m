%无状态约束优化问题测试
warning off
clear 
clc
global Lt La1 La T;
global f fx fu L Lx Lu cL cLx cLu n lx lu Nh NO;  

NO=3;%NO+1为有限单元上取的配置点数,配置点适当增加，则求解越精确，但是求解时间越长
Nh=0.1;%有限元长度，有限单元数过大会导致微分方程模拟失败
clc
belta=0.33;gama=0.18;k=0.19;
 f=@(x,u)[-(1-u(1))*belta*x(1)*x(3),...
         (1-u(1))*belta*x(1)*x(3)-k*x(2),...
          k*x(2)-gama*x(3)];
 fx=@(x,u)[-(1-u(1))*belta*x(3),0,-(1-u(1))*belta*x(1);...
           (1-u(1))*belta*x(3),-k,(1-u(1))*belta*x(1);...
           0,k,-gama];
 fu=@(x,u)[belta*x(1)*x(3);-belta*x(1)*x(3);0];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 L=@(x,u)(1-u(1))*belta*x(1)*x(3);
 Lx=@(x,u)[(1-u(1))*belta*x(3),0,(1-u(1))*belta*x(1)];
 Lu=@(x,u)-belta*x(1)*x(3);
%%%%%%%%%%%%%%%%%
x0=[1000,100,10];%可适当取值
x0=10*x0/sum(x0);
n=10;
lu=1;
%%%%%%%%%%%%%%%%%%%%%对控制的设置
Ul=zeros(n,lu);
UM=0.33*ones(n,lu);
U=UM/2+Ul/2;
% U=UM;


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
