function [J,dJu,XX,XU]=funocfe(U,x0);
global Lt La1 La T;
global f fx fu L Lx Lu cL cLx cLu n lx lu Nh NO
U=reshape(U,lu,n)';
m=NO;h=Nh;%h为单元长度，m+1为有限单元上的配置点数目,总长tf=h*n.
%本程序中进行均匀插值，即T=h:h/(m+1):h;
%[Lt,La1,La,T]=OCFElagrange(linspace(0,h,m+2));%进行求解时需要的拉格朗日插值函数,La1为t=1时
lis=Lt;lis(:,1)=[];lis(1,:)=[];lis1=inv(lis');bt=lis1(end,:);

%x0=[1000,1,1];
lx=length(x0);%定义状态初值
sizeX1=m+1;%定义optfun时要用到的函数

XX=[];
Flag=[];ode=[];II=[];%过程值存储
li=(m+1)*lx*n;
ODEX=sparse(li,li);ODEU=sparse(li,n*lu);
J=0;Jx=sparse(1,li);Ju=sparse(1,n*lu);
dJu=sparse(1,n*lu);

for k=1:n;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u=U(k,:);fi=@(x)f(x,u);fix=@(x)fx(x,u);fiu=@(x)fu(x,u);
optfun=@(X)OCFEode4(X,x0,2,sizeX1,fi,fix,fiu,lu);%每一单元上的微分约束及其雅克比矩阵,由于Lt已经规范化，此处取h=2
li=4*m+4;dh=h/li;X=x0;lix=x0;
for i=1:li;lix=lix+dh*f(lix,u);X=[X;lix];end;X=interp1(0:dh:h,X,T);%猜测一个初始值方便用牛顿迭代法求解
%上一步用欧拉迭代法结合插值对有限元上插值点状态变量进行预测，求得一个合适的迭代初始点
%这一步是十分必要的，否则，不但求解时间将增加到原来的三倍，更重要的是解是错误的！！！！！
X=reshape(X',numel(X),1);%将原来的（m+1)*lx矩阵转换为列向量形式，以符合@optfun函数的输入规范
flag=0;
for i=1:20;%牛顿迭代法求解有限元上的微分约束方程，得到m+1个插值点
 %[li,A]=optfun(X);
if i==1;[li,A]=optfun(X);else li=optfun(X);end;%使用定点牛顿迭代方法减少计算雅克比矩阵的时间
if norm(li)<=1e-6;flag=1;break;end;
dX=-A\li;X=X+dX;
end;
Flag=[Flag,flag];ode=[ode,norm(li)];II=[II,i];
Xr=reshape(X,lx,sizeX1);Xr=Xr';%将列向量X还原为（m+1)*lx矩阵形式Xr
XX=[XX;Xr];%储存状态变量
x0=Xr(end,:);

%计算ODEX，ODEU
if nargout>1;
[li,odex,odex0,odeu]=optfun(X);
li0=(m+1)*(k-1)*lx;li1=li0+1:li0+(m+1)*lx;li2=li0+1:li0+(m+1)*lx;
li3=li0-lx+1:li0;li4=(k-1)*lu+1:(k-1)*lu+lu;
ODEX(li1,li2)=odex;
if k>1;ODEX(li1,li3)=odex0;end;
ODEU(li1,li4)=odeu;
end;
%计算目标函数
li0=(m+1)*(k-1)*lx;
Ji=0;liju=sparse(1,lu);
for i=1:m+1;
    xi=Xr(i,:);%Xr, m*lx
    bi=bt(i);
    Ji=Ji+bi*L(xi,u);
    if nargout>1;Jx(li0+(i-1)*lx+1:li0+i*lx)=bi*Lx(xi,u);
    liju=liju+bi*Lu(xi,u);end;
end;
if nargout>1;Ju((k-1)*lu+1:k*lu)=liju;end;
J=J+Ji;

end;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout>1;
XU=-ODEX\ODEU;
dJu=Jx*XU+Ju;
dJu=dJu';%输出列向量
dJu=full(dJu);
end;

%Uli=U;li=20;ep=1e-5;Uli(li)=U(li)+ep;[li0,sb]=funOCFE(U);li1=funOCFE(Uli);
%[sb(li),(li1-li0)/ep]
