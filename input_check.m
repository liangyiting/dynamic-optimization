%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%对所提供的输入变量进行检查，主要是检查梯度计算是否准确,以及所给的形式是否正确%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check=0;
ep=1e-4;
u=1/2*Ul(1,:)+1/2*UM(1,:);
lx=length(x0);lu=length(u);
if size(x0,1)~=1;x0,disp('初始值必须以行向量的形式输入');check=1;end;
if size(f(x0,u),2)~=lx|size(f(x0,u),1)~=1;f(x0,u),disp('函数f的输出必须是(1,lx)矩阵');check=1;end;
if size(fx(x0,u),2)~=lx|size(fx(x0,u),1)~=lx;f(x0,u),disp('函数fx的输出必须是(lx,lx)矩阵');check=1;end;
if size(fu(x0,u),2)~=lu|size(fu(x0,u),1)~=lx;f(x0,u),disp('函数fu的输出必须是(lx,lu)矩阵');check=1;end;
if size(Lx(x0,u),2)~=lx|size(Lx(x0,u),1)~=1;Lx(x0,u),disp('函数Lx的输出必须是(1,lx)矩阵');check=1;end;
if size(Lu(x0,u),2)~=lu|size(Lu(x0,u),1)~=1;Lu(x0,u),disp('函数Lu的输出必须是(1,lu)矩阵');check=1;end;

f0=f(x0,u);L0=L(x0,u);
df0=fx(x0,u);dL0=Lx(x0,u);
try c0=cL(x0,u);dc0=cLx(x0,u);dc=[];end;
dL=[];df=[];
for i=1:length(x0);
    xi=x0;xi(i)=x0(i)+ep;
    fi=f(xi,u);Li=L(xi,u);
    try ci=cL(xi,u);dc(:,i)=((ci-c0)/ep)';end;
    df(:,i)=([fi-f0]/ep)';
    dL(:,i)=([Li-L0]/ep)';
end;
if norm(df-df0)>1e-9*norm(df);'f或者fx函数输入有误',check=1;end
if norm(dL-dL0)>1e-9*norm(dL);'L或者Lx函数输入有误',check=2;end
try if (norm(dc-dc0))>1e-9*norm(dc);'cL或者cLx函数输入有误',check=3;end;end;
    
ep=1e-4;
u=1/2*Ul(1,:)+1/2*UM(1,:);
f0=f(x0,u);L0=L(x0,u);
df0=fu(x0,u);dL0=Lu(x0,u);
try c0=cL(x0,u);dc0=cLu(x0,u);dc=[];end;
dL=[];df=[];
for i=1:length(u);
    ui=u;ui(i)=u(i)+ep;
    fi=f(x0,ui);Li=L(x0,ui);
    try ci=cL(x0,ui);dc(:,i)=((ci-c0)/ep)';end;
    df(:,i)=([fi-f0]/ep)';
    dL(:,i)=([Li-L0]/ep)';
end;
if norm(df-df0)>1e-9*norm(df);'f或者fx函数输入有误',check=4;end
if norm(dL-dL0)>1e-9*norm(dL);'L或者Lx函数输入有误',check=5;end
try if (norm(dc-dc0))>1e-9*norm(dc);'cL或者cLx函数输入有误',check=6;end;end;
if lu~=size(U,2)|n~=size(U,1);'U输入维数不对，U应该是（n，lu）维矩阵',check=7;end;
if lu~=size(Ul,2)|n~=size(Ul,1);'Ul输入维数不对，Ul应该是（n，lu）维矩阵',check=8;end;
if lu~=size(UM,2)|n~=size(UM,1);'UM输入维数不对，UM应该是（n，lu）维矩阵',check=9;end;
if sum(sum(Ul>=UM))>0;'下界Ul中某一位置上的值大于上界UM中同一位置上的值，请检查',check=10;end;
if check>0;'检查没有通过，请按照提示检查输入',error,end;