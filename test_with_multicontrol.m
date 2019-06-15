%多控制变量优化问题测试
%%%%%%%%%%%%%%%%问题描述
%min J=int_0^T((1-u(1))*beta*S*I)dt
%dS/dt=-(1-u(1))*beta*S*I-u(2)*S,
%dE/dt=(1-u(1))*beta*S*I-muE*E,
%dI/dt=muE*E-muI*I-u(3)*I
warning off
clear 
clc
global Lt La1 La T;
global f fx fu L Lx Lu cL cLx cLu n lx lu Nh NO;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NO=1;%NO+1为有限单元上取的配置点数,配置点适当增加，则求解越精确，但是求解时间越长
Nh=1;%有限元长度，有限单元数过大会导致微分方程模拟失败，应当选<20为好
clc

%%%%%%%%%%%%%%参数设置
beta=0.33;muE=0.19;muI=0.18;
f=@(x,u)[-(1-u(1))*beta*x(1)*x(3)-u(2)*x(1),...
         (1-u(1))*beta*x(1)*x(3)-muE*x(2),...
        muE*x(2)-muI*x(3)-u(3)*x(3)
        ];
fx=@(x,u)[-(1-u(1))*beta*x(3)-u(2),0,-(1-u(1))*beta*x(1);
          (1-u(1))*beta*x(3),-muE,(1-u(1))*beta*x(1);
          0,muE,-muI-u(3)
        ];
fu=@(x,u)[beta*x(1)*x(3),-x(1),0;-beta*x(1)*x(3),0,0;0,0,-x(3)];
%%%%%%%%%%%%%%%%%%%%%%%目标函数L(x,u)
L=@(x,u)((1-u(1))*beta*x(1)*x(3));
Lx=@(x,u)[(1-u(1))*beta*x(3),0,(1-u(1))*beta*x(1)];
Lu=@(x,u)[-beta*x(1)*x(3),0,0];
%%%%%%%%%%%%%%%%%%%%%%模拟时间分段数n，时间段长度h.  注：对不同的tf，可调整n的值以减少计算时间
x0=[1 .1 .01];
n=100;%模拟长度
lu=3;
%%%%%%%%%%%%%%%%%%%%%对控制的设置
Ul=zeros(n,lu);
ua=beta;ub=1;uc=1;UM=[ua*ones(n,1),ub*ones(n,1),uc*ones(n,1)];%3双控制；控制变量下界um和上界uM
U=UM/2+Ul/2;


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




