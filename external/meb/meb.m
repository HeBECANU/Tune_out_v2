function [c r] = meb(X,mode)
% This is a function for calculate the smallest enclosing circle for a
% set of points. User can choose a way to do the job.
%==========================================================================
% Mode:
% 'cvx':    Calculate the circle by transforming the problem to its convex form
% 'elzinga':Calculate the circle by using Elzinga-Hearn algorithme
% 'emo':    Calculate the circle by using Emo Welzl algorithme.
%
% X : 
% The set of points for calculation. The size should be N*2(N points)
% You might need to change the maximum recursive time to ensure the functionality,
% Command matlab: set(0,'RecursionLimit',1000);
% /*1000 will be sufficient for a thousand points.*/
% 
% Attention:
% With the mode 'elzinga', you can only calculate for less than
% 50 points, otherwise the matlab will be turned down.
% 
% c : The certre of output circle with size(c)=[1,2]
% r : The radius of output circle
%==========================================================================
% An exemple:
%   clear all
%   n=10;
%   dim=2;
%   X = randn(n,dim) + 1.5*ones(n,1)*[1 2.5];
%   [c r] = meb(X,'emo');
%   ang=0:0.01:2*pi; 
%   xp=r*cos(ang);
%   yp=r*sin(ang);
%   figure
%   plot(X(:,1),X(:,2),'r*');
%   hold on
%   plot(c(1)+xp,c(2)+yp);
%
%
%  Yi Yu
%  yi.yu@insa-rouen.fr
%  INSA Rouen, France
%  12, 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(mode,'cvx'))
    [c,r]=cvx_MEB(X);
else 
    if(strcmp(mode,'elzinga'))
        [c,r]=elzinga_MEB(X);
    else
        if(strcmp(mode,'emo'))
            [c,r]=emo_MEB(X);
        else
            display('This mode does not exist, help meb for more information');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cvx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,r]=cvx_MEB(X)
% Version QP of problem : min R^2 with || x_i - c||^2 <= R^2 ;  i=1,n 
nb=size(X,1);
M=X*X';
Nx=diag(M);
cvx_begin quiet
cvx_begin
   cvx_precision best;
   variable x(nb) nonnegative;
   dual variable rho;
   minimize(x'*M*x-Nx'*x);
   subject to
        rho:sum(x)==1;
cvx_end

c=X'*x;
r=sqrt(x'*M*x-rho);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% elzinga %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,r]=elzinga_MEB(X)
    p1=X(1,:);
    p2=X(2,:);
    [c,r]=elzinga_etape2(p1,p2,X);
end

function [c,r]=elzinga_etape2(p1,p2,X)
    n=size(X,1);
    [c,r]=cercle2points(p1,p2);
    d=diag((X-ones(n,1)*c)*(X-ones(n,1)*c)');
    if(max(d)> r^2)
        pk=X(d==max(d),:);
        [c,r]=elzinga_etape3(p1,p2,pk,X);
    end
end

function [c,r]=elzinga_etape3(p1,p2,pk,X)
    n=size(X,1);
    if (typeTriangle(p1,p2,pk)<=0)
            [c1,r1]=cercle2points(p2,pk);
            [c2,r2]=cercle2points(p1,pk);
            if(r1>=r2)
                [c,r]=elzinga_etape2(p2,pk,X);
            else
                [c,r]=elzinga_etape2(p1,pk,X);
            end
    else        
        [c,r]=cercleTriangle(p1,p2,pk);
        d=diag((X-ones(n,1)*c)*(X-ones(n,1)*c)');
        if(max(d)>r^2)
             [c,r]=elzinga_etape4(c,p1,p2,pk,d,X);
        end
    end
end

function [c,r]=elzinga_etape4(c,p1,p2,pk,d,X)
    pl=X(d==max(d),:);
    pl=pl(1,:);
    pn=[p1;p2;pk];
    dpl=diag((pn-ones(3,1)*pl)*(pn-ones(3,1)*pl)');
    q=pn(dpl==max(dpl),:);
    ps=[pn;pl];
    %Verifier la position des points selon le droit qui passe c et q
    k=(c(2)-q(2))/(c(1)-q(1));   
    b=c(2)-k*c(1);
    position = sign(k.*(ps(:,1))+b-ps(:,2));
    for i=1:3
       if(sign(position(4))*sign(position(i))==-1)
              p1=q;
              p2=ps(i,:);
              pk=pl;
       end
    end
    [c,r]=elzinga_etape3(p1,p2,pk,X);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% emo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,r]=emo_MEB(X)
nb=size(X,1);
%B: points set contains those on the circle
B=[];
%C: centre du cercle
C=[];
b=0;
[c,r,C]=emo_minidisk(X,nb,B,b,C);
end
    
function [c,r,C]=emo_minidisk(P,p,B,b,C)
%P: points set contains those inside the circle
%p: indice of the points picked(always the last row of P)
%B: points set contains those on the circle
%b: card(B),0<=b<=3
if(p==0 || b==3)
     switch b
          case 0
            c=[];
            r=0;
          case 1
            c=B;
            r=0;
          case 2
            [c,r]=cercle2points(B(1,:),B(2,:));
            C=[C;c];
          case 3
           [c,r]=cercleTriangle(B(1,:),B(2,:),B(3,:));
            C=[C;c];
     end
     return 
end
[c,r,C]=emo_minidisk(P,p-1,B,b,C);
C=[C;c];

if isempty(c) || (sqrt((P(p,1)-c(1)).^2+(P(p,2)-c(2)).^2)>r)
   B=[B;P(p,:)];
   [c,r,C]=emo_minidisk(P,p-1,B,b+1,C);
   C=[C;c];
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%% other functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate circle from 2 points
function [c,r]=cercle2points(p1,p2)
c=(p1+p2)/2;
r=sqrt((p1(1)-c(1))^2+(p2(2)-c(2))^2)/2;
end

% Calculate circle from 3 points
function [c,r]=cercleTriangle(p1,p2,p3)
cen1=(p1+p2)/2;        
cen2=(p2+p3)/2;         

k1=-1/((p1(2)-p2(2))/(p1(1)-p2(1)));   
b1=cen1(2)-k1*cen1(1);

k2=-1/((p2(2)-p3(2))/(p2(1)-p3(1)));    
b2=cen2(2)-k2*cen2(1);

x0=-(b1-b2)/(k1-k2);             
y0=-(-b2*k1+b1*k2)/(k1-k2);
c=[x0 y0];                                    
r=sqrt((y0-p1(2))^2+(x0-p1(1))^2); 
end        
end

