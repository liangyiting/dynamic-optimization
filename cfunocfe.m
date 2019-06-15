function [ck,ceq,dck,dceq]=cfunocfe(U,x0)
global f fx fu L Lx Lu cL cLx cLu n lx lu Nh NO
[~,~,X1,XU]=funocfe(U,x0);
U=reshape(U,lu,n)';
m=size(X1,1)/n-1;
ck=[];dck=sparse(n*(m+1),size(XU,2));%dck(i,j)=ci/uj;而fmincon函数中需要的是dck(i,j)=cj/ui
for k=1:n*(m+1);
    x=X1(k,:);li0=ceil(k/(m+1));u=U(li0,:);
    li=cL(x,u);ck=[ck;li];
    if nargout>1;
    dlix=cLx(x,u);dliu=cLu(x,u);
    li2=sparse(1,n*lu);li2((li0-1)*lu+1:(li0-1)*lu+lu)=dliu;
    li=lx*(k-1);lixu=XU(li+1:li+lx,:);
    lidck=dlix*lixu;lidck=lidck+li2;
    dck(k,:)=lidck;
    end;
end;
ceq=[];
dceq=[];
dck=dck';
