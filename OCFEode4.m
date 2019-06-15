function [a,b,odex0,odeu]=OCFEode4(X,x0,h,m,fi,fix,fiu,lu);

global Lt;
if nargout==1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X=reshape(X,length(x0),m);
X=X';   
odex=[];
lx=length(x0);
for k=1:m;
    x=X(k,:);
    xt=2/h*Lt(:,k+1)'*[x0;X];
    ode(k,:)=xt-fi(x);
end;
a=ode;a=a';
a=reshape(a,numel(a),1);
elseif nargout==2;%%%%%%%%%%%%%%%%%%%%%%%%%
X=reshape(X,length(x0),m);
X=X';   
odex=[];
lx=length(x0);
for k=1:m;
    x=X(k,:);
    xt=2/h*Lt(:,k+1)'*[x0;X];
    ode(k,:)=xt-fi(x);
    liodex=[];
    for i=1:m;
    if i==k;odexi=Lt(i+1,k+1)*eye(lx)*2/h-fix(x);else odexi=Lt(i+1,k+1)*eye(lx)*2/h;end;
    liodex=[liodex,odexi];
    end;
    odex=[odex;liodex];
end;
a=ode;b=odex;
a=a';
a=reshape(a,numel(a),1);
elseif nargout==4; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X=reshape(X,length(x0),m);
X=X';   
odex=[];
lx=length(x0);
for k=1:m;
    x=X(k,:);
    xt=2/h*Lt(:,k+1)'*[x0;X];
    ode(k,:)=xt-fi(x);
    liodex=[];
    for i=1:m;
    if i==k;odexi=Lt(i+1,k+1)*eye(lx)*2/h-fix(x);else odexi=Lt(i+1,k+1)*eye(lx)*2/h;end;
    liodex=[liodex,odexi];
    end;
    odex=[odex;liodex];
    li1=(k-1)*lx;odex0(li1+1:li1+lx,1:lx)=2/h*Lt(1,k+1)*eye(lx);
    odeu(li1+1:li1+lx,1:lu)=-fiu(x);
end;
a=ode;b=odex;
a=a';
a=reshape(a,numel(a),1);
end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
