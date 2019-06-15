%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% 输入变量通过检查后，进行下面的计算%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lx=length(x0);
sb=[];
m=NO;
h=Nh;
tf=n*h;%末端时间
[Lt,La1,La,T]=OCFElagrange(linspace(0,h,m+2));%进行求解时需要的拉格朗日插值函数,La1为t=1时
pro=@(u)reshape(u,lu,length(u)/lu)';
ipro=@(u)reshape(u',numel(u),1);
U=ipro(U);Ul=ipro(Ul);UM=ipro(UM);
fun=@(U)funocfe(U,x0);[J,dJu,XX,XU]=fun(U); 
L(x0,U(:,1));
try cL(x0,U(:,1));%如果定义了逐点状态不等式约束
    cfun=@(U)cfunocfe(U,x0);[ck,sssss,dck,ssss]=cfun(U);%梯度矩阵格式 dck(i,j)=cj/ui
catch cfun=[];%如果没有定义逐点状态不等式约束
end;

%求解设置，可选的算法有内点法'interior-point','sqp','trust-region-reflective' (current default),
%'active-set''interior-point' (will become default in a future release),'sqp'
%优化过程显示
% @optimplotx plots the current point
% @optimplotfunccount plots the function count
% @optimplotfval plots the function value
% @optimplotconstrviolation plots the maximum constraint violation
% @optimplotstepsize plots the step size
% @optimplotfirstorderopt plots the first-order optimality measure
options=optimset('display','iter','MaxFunEvals',maxfunevals,'GradConstr','on','GradObj','on',...
    'Algorithm',algorithm,'PlotFcns',plotfcns);
[Uopt,fval,exitflag,output,lambda]= fmincon(fun,U,[],[],[],[],Ul,UM,cfun,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[J,dJu,X1,XU]=fun(Uopt);
Xopt=[x0;X1];ttx=0:h/(m+1):tf;
Uopt=pro(Uopt);ttu=0:h:tf-h;
figure(2);plot(ttu,Uopt);title('控制轨迹');xlabel('时间');
figure(3);plot(ttx,Xopt);title('状态轨迹');xlabel('时间');
