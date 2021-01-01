function [x,ccode,it,er] = linppsr(A,b,numit,epsilon,ncomp)
%   LINPPSR - optimal (exact) componentwise estimation of the united solution 
%   set to interval linear system of equations.
%  
%   x = linppsr(A,b) computes optimal componentwise lower and upper estimates  
%   of the solution set to interval linear system of equations Ax = b,
%   where A - square interval matrix, b - interval right-hand side vector. 
%   
%  
%   [x,ccode,it,er] = linppsr(A,b,numit,epsilon,ncomp) computes vector x of 
%   optimal componentwise estimates of the solution set to interval linear 
%   system Ax = b with accuracy no more than epsilon and after the number of   
%   iterations no more than numit. Optional input argument ncomp indicates 
%   a component's number of interval solution in case of computing the estimates 
%   for this component only. If this argument is omitted, all componentwise 
%   estimates is computed. 
%   Additional output arguments contain information:
%        ccode - termination code: 
%                0 - optimal (exact) estimation of solution set is obtained,
%                1 - estimates of solution set is computed with prescribed 
%                    accuracy epsilon;
%                2 - estimates of solution set is computed in maximum 
%                    number of iterations numit; 
%        it    - number of really performed iterations;
%        er    - accuracy of computed estimates.
   
%   Procedure LINPPSR uses adaptive Partitioning of the Parameter Set 
%   (PPS algorithm) described in
%
%   S.P. Shary
%   Parameter partition methods for optimal numerical solution of interval 
%   linear systems, in: Computational Science and High-Performance Computing III. 
%   The 3rd Russian-German advanced research workshop, Novosibirsk, Russia, 
%   July 23-27, 2007, E. Krause, Yu.I. Shokin, M. Resch, N.Yu. Shokina, eds. 
%   - Berlin-Heidelberg, Springer, 2008. - P. 184-205.
%  
%   Example 1. 
%   intval A = 
%   [    2.0000,    4.0000] [   -2.0000,    1.0000] 
%   [   -1.0000,    2.0000] [    2.0000,    4.0000] 
%   intval b = 
%   [   -2.0000,    2.0000] 
%   [   -2.0000,    2.0000] 
%  
%   [x,ccode,it] = linppsr(A,b)
%   intval x = 
%   [   -4.0000,    4.0000] 
%   [   -4.0000,    4.0000] 
%   ccode =
%           0     0
%           0     0
%   it =
%           2     2
%           2     2
%   Optimal (exact) estimation of solution set is obtained (ccode = 0), 
%   every componentwise estimate is computed in 2 iterations (it = 2).
%  
%   Example 2.
%   intval A = 
%   [   10.0000,   10.0000] [    0.0000,    2.0000] [    0.0000,    2.0000] 
%   [    0.0000,    2.0000] [   10.0000,   10.0000] [    0.0000,    2.0000] 
%   [    0.0000,    2.0000] [    0.0000,    2.0000] [   10.0000,   10.0000] 
%   intval b = 
%   [   -1.0000,    1.0000] 
%   [   -1.0000,    1.0000] 
%   [   -1.0000,    1.0000] 
%  
%   [x,ccode,it,er] = linppsr(A,b,realmax,0.02,1)
%   intval x = 
%   [   -0.1524,    0.1524] 
%   ccode =
%          1     1
%   it =
%          2     2
%   er =
%         0.0081    0.0081
%   Lower and upper estimates of the first component of interval solution (ncomp = 1) 
%   is computed with  prescribed accuracy epsilon = 0.02.
% 
%   Example 3.
%   intval A = 
%   [    3.0,    5.0] [   -0.6,    0.4] [   -0.6,    0.4] [   -0.6,    0.4] 
%   [   -0.6,    0.4] [    3.0,    5.0] [   -0.6,    0.4] [   -0.6,    0.4] 
%   [   -0.6,    0.4] [   -0.6,    0.4] [    3.0,    5.0] [   -0.6,    0.4] 
%   [   -0.6,    0.4] [   -0.6,    0.4] [   -0.6,    0.4] [    3.0,    5.0]  
%   intval b = 
%   [   -3.0,    3.0] 
%   [   -3.0,    3.0] 
%   [   -3.0,    3.0] 
%   [   -3.0,    3.0] 
%
%   [x,ccode,it] = linppsr(A,b,5,0,1)
%   intval x = 
%   [   -2.5001,    2.5001] 
%   ccode =
%           2     2
%   it =
%           5     5
%   Lower and upper estimates of the first component of interval solution (ncomp = 1)  
%   in 5 iterations (it = 5).
%
%   [x,ccode,it] = linppsr(A,b)
%   intval x = 
%   [   -2.5000,    2.5000] 
%   [   -2.5000,    2.5000] 
%   [   -2.5000,    2.5000] 
%   [   -2.5000,    2.5000] 
%   ccode =
%           0     0
%           0     0
%           0     0
%           0     0
%   it =
%           9     9
%           9     9
%           9     9
%           9     9
%  
%   Optimal (exact) estimation of solution set is obtained (ccode = 0), 
%   every componentwise estimate is computed in 9 iterations (it = 9). 
%   
%   Remark. Run time of the algorithm depends on the dimension of the interval 
%   linear algebraic system (ILAS) to be solved and properties of its matrix. 
%   If the system has a large dimension or its matrix is near the bounds of 
%   singularity, then the run time may be very large too. To save time, you 
%   may wish to compute lower and upper estimates of only one component of 
%   the solution set and allot a certain number of iterations for doing that. 
%   This can be set by the parameter 'it'. 
%   In case of acceptable run time, optimal estimates of the solution set 
%   can be obtained by the use of procedure [x,ccode,it,er] = linppsr(A,b).
%   Otherwise, one computes outer componentwise estimates with the accuracy
%   'er' in the prescribed number of iterations, and the estimates is not 
%   necessarily optimal. 
%   
%   Dmitry Lyudvin, Sergey P. Shary, 2010.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   checking correctness of the data inputted 
  
[n,n1] = size(A);
m1 = length(b);
  
if n~=n1
    error('Matrix is not square')
end
if n~=m1
    error('Inconsistent dimensions of matrix and right-hand side vector')
end
  
A = intval(A);
b = intval(b);

if det(mid(A))==0
    error('Matrix is singular')
end
  
if max(abs(eig(abs(inv(mid(A)))*rad(A))))>=1
    error('Matrix is singular')
end  
InvA = inv_HBR(A);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
szL = 3;                         %   szL - length of the list of active records
rt = 0.05;
empty = infsup(NaN,NaN);
dQ = repmat(empty,n,n);
sv = zeros(2,n);
  
if nargin==2
    epsilon = 0;
    numit = realmax;
end
if (nargin==3)
    epsilon = 0;
end
if (nargin==5)
    nmin = ncomp;
    nmax = ncomp;
    x = repmat(empty,1,1);
    ccode = zeros(1,2);
    it = zeros(1,2);
    er = zeros(1,2);
else
    nmin = 1;
    nmax = n;
    x = repmat(empty,n,1);
    ccode = zeros(n,2);
    it = zeros(n,2);
    er = zeros(n,2);
end
if epsilon <= 0 && numit == realmax
    mode = 1;
else mode = 2;
end
  
for m = 1:2
    if m==2                         %   if m = 1 then compute the lower estimates
        b = -b;                     %   if m = 2 then compute the upper estimates
    end
  
    for nu = nmin:nmax
        vx = HBR(A,b);              %  outer interval estimates of solution set 
        if all(isnan(vx))
            error('The system is incompatible');
        end
        omega = Inf; 
                                    %   first (leading) record in the list:
        L(1).gamma = inf(vx(nu));   %   lower endpoint of nu-th component
                                    %   of interval estimate
        L(1).Q = A;                 %   interval matrix of the system
        L(1).r = b;                 %   interval right-hand side vector
        L(1).Y = InvA;              %   inverse interval matrix
        L(1).x = vx;                %   outer interval estimates of solution set
             
        %   initial (zeros) check matrix and their recording in list:
        L(1).W = zeros(n,n);        %   matrix W 
        L(1).s = zeros(n,1);        %   column vector s
        L(1).t = zeros(1,n);        %   row vector t
        
        Lk = [];
        cstL = Inf;
        Nit = 0;
        if Nit == 0
            cont = 1;
        else
            cont = any(any(rad(L(1).Q))~=0);
        end
        while cont && abs(omega-L(1).gamma)>epsilon && Nit<=numit
            Omega = [];
            K = [];
            Wch = false;
            sch = false;
  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
            %   monotony test
            for i = 1:n
                dQ(i,1:n) = -L(1).Y(nu,i)*L(1).x(1:n);  %   dQ - enclosures of derivatives
            end                             %   of the solution with respect to the matrix 
  
            for i = 1:n
                for j = 1:n
                    if dQ(i,j)>=0 
                        L(1).Q(i,j) = inf(L(1).Q(i,j)); 
                        L(1).W(i,j) = 1;
                    end
                    if dQ(i,j)<=0 
                        L(1).Q(i,j) = sup(L(1).Q(i,j));  
                        L(1).W(i,j) = -1;
                    end
                end
            end
  
            for i = 1:n             %   dR = L(1).Y(nu,:) - interval enclosures of derivatives
                                    %   of the solution with respect to the right-hand side vector                                    
                if L(1).Y(nu,i)>=0 
                    L(1).r(i) = inf(L(1).r(i)); 
                    L(1).s(i) = -1;
                end
                if L(1).Y(nu,i)<=0  
                    L(1).r(i) = sup(L(1).r(i));  
                    L(1).s(i) = 1;
                end
            end
  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
            %   find the element providing the maximal product of the width 
            %                                  by the derivative estimate
            %   z1 - computed element of the matrix
            %   im1 jm1 - its indexes
            [c,I] = max(mag(dQ).*diam(L(1).Q));
            [z1,jm1] = max(c);
            im1 = I(jm1);  
  
            %   z2 - computed element of the right-hand side vector
            %   im2 - its index
            [z2,im2] = max(mag(L(1).Y(nu,1:n)).*(diam(L(1).r))'); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     refining W by searching 2õ2-submatrix 
            out = false;
            for i = 1:n 
              for j = 1:n
                if (i~=im1 && j~=jm1) 
                  if (L(1).W(i,j)~=0 && L(1).W(i,jm1)~=0 && L(1).W(im1,j)~=0)
                    L(1).W(im1,jm1) = L(1).W(i,j)*L(1).W(i,jm1)*L(1).W(im1,j);
                    out = true;
                    break
                  end
                end
              end
              if out 
                break 
              end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   generating the systems-descendants   
            W = L(1).W;
            s = L(1).s;
            t = L(1).t;
            Q = L(1).Q;
            r = L(1).r;
  
            if z1 > z2
                switch W(im1,jm1)
                case 0     %  case of two systems-descendants Q' è Q"
                    Q(im1,jm1) = inf(Q(im1,jm1));
                    W(im1,jm1) = 1;
                    Omega{length(Omega)+1} = [im1 jm1];
                    Wch = true;
                    m1 = split(1,1);
                    W = L(1).W;
                    s = L(1).s;
                    t = L(1).t;
                    Omega = [];
                    K = [];
  
                    Q(im1,jm1) = sup(L(1).Q(im1,jm1));
                    W(im1,jm1) = -1;
                    Omega{length(Omega)+1} = [im1 jm1];
                    Wch = true;
                    m2 = split(1,1);
                    mu = min(m1(nu),m2(nu));
                case 1     %   case of one system-descendant Q'  
                    Q(im1,jm1) = inf(Q(im1,jm1));
                    m1 = split(0,1);
                    mu = m1(nu);
                case -1    %  case of one system-descendant Q"
                    Q(im1,jm1) = sup(Q(im1,jm1));
                    m1 = split(0,1);
                    mu = m1(nu);
                end
            else
                switch s(im2)
                case 0  %  case of two systems-descendants r' è r"  
                    r(im2) = inf(r(im2));
                    s(im2) = -1;
                    K(length(K)+1) = im2;
                    sch = true;
                    m1 = split(1,0);
                    s = L(1).s;
                    W = L(1).W;
                    t = L(1).t;
                    Omega = [];
                    K = [];
  
                    r(im2) = sup(L(1).r(im2));
                    s(im2) = 1;
                    K(length(K)+1) = im2;
                    sch = true;
                    m2 = split(1,0);
                    mu = min(m1(nu),m2(nu));
                case 1     %  case of one system-descendant r'
                    r(im2) = sup(r(im2));
                    m1 = split(0,0);
                    mu = m1(nu);
                case -1    %  case of one system-descendant r" 
                    r(im2) = inf(r(im2));
                    m1 = split(0,0);
                    mu = m1(nu);
                end 
            end
  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
            L(1) = [];              %   delete the leading record from the list
            if  omega>mu  
                omega = mu;         %   compute the new value of omega
            end
                                    %   delete unpromising records
            tL = length(L);         
            tLk = length(Lk);       
  
            if mode==1 && ~isempty(L) && ~any(any(rad(L(1).Q)))
                for i = 1:tL
                    Lk = transfer(Lk,L,tLk+i,i);
                end
                L = [];
            end
  
            if isempty(L) && ~isempty(Lk)       
                i = 1;                          
                while i<=length(Lk)
                    if Lk(i).gamma>omega
                        Lk(i) = [];
                    else i = i+1;
                    end
                end
                tLk = length(Lk);
                if tLk~=0
                    [H,I] = sort([Lk.gamma]);   
                    if mode == 1
                        T = min(szL,tLk); 
                        for j = 1:T
                            L = transfer(L,Lk,j,I(j));
                        end
                    else
                        cstL = (rt*omega+Lk(I(1)).gamma)/(1+rt); 
                        T = 1;
                        while H(T)<cstL
                            L = transfer(L,Lk,T,I(T));
                            T = T+1;
                            if T>tLk
                                break
                            end
                        end
                        T=T-1;
                    end
                    J = sort(I(1:T));
                    for j = 1:T
                        Lk(J(j)-j+1) = [];
                    end
                end
            end
            tL = length(L);
            if ~isempty(L) && isempty(Lk) 
                for i = 1:tL
                    if L(i).gamma>omega
                        L(i:tL) = []; 
                        break
                    end
                end
            end
            Nit = Nit+1;
            if isempty(L)&& isempty(Lk)
                ocenka = omega;
                break
            end
            ocenka = L(1).gamma;
            cont = any(any(rad(L(1).Q))~=0);
        end
        sv(m,nu) = ocenka;   %  estimate of the nu-th component of solution set 
        if (nargin==5)
            e = 1;
        else
            e = nu; 
        end
        it(e,m) = Nit-1;
        er(e,m) = abs(omega-ocenka);
        if er(e,m)<=epsilon
            ccode(e,m) = 1;
            if er(e,m)==0
                ccode(e,m) = 0;
            end  
        end
        if Nit>=numit
            ccode(e,m) = 2;
        end
  
        Lk = []; 
        L = []; 
    end
end
  
if (nargin==5)
    x = infsup(sv(1,nu),-sv(2,nu));         %   estimate of the nu-th component 
else
    for nu = nmin:nmax
        x(nu) = infsup(sv(1,nu),-sv(2,nu)); %   interval enclosure 
    end 
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function Mn = split(ind,fl)  
%   procedure of generating the systems-descendants 
invQ = L(1).Y;
vx = HBR(Q, r);
gamma = inf(vx(nu));
Mn = mid(Q)\mid(r);
  
%   if gamma < omega and fl = true then recalculate W, s, t   
if gamma<omega && ind
    tch = false;
    Wc = false;
    sc = false;
    tc = false;
    Lambda = [];
    while (Wch || sch || tch)
        if (sch || tch)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   recalculation of W
            if sch
                for i = 1:length(K)
                    for j = 1:n
                        if t(j)~=0
                            a = s(K(i))*t(j);
                            if W(K(i),j)~=a
                                Wc = true;
                                Omega{length(Omega)+1} = [K(i) j];
                            end  
                            W(K(i),j) = a;
                        end
                    end
                end
            end
  
            if tch
                for i = 1:length(Lambda)
                    for j = 1:n
                        if s(j)~=0
                            a = s(j)*t(Lambda(i));
                                if W(j,Lambda(i))~=a
                                    Wc = true;
                                    Omega{length(Omega)+1} = [j Lambda(i)];
                                end 
                                W(j,Lambda(i)) = a;
                        end
                    end
                end
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        if (Wch || tch)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   recalculation of s  
            if Wch
                for i = 1:length(Omega)
                    i1 = Omega{i}(1);
                    i2 = Omega{i}(2);
                    if t(i2)~=0
                        a = W(i1,i2)/t(i2);
                        if s(i1)~=a
                            sc = true;
                            K(length(K)+1) = i1;
                        end 
                        s(i1) = a;
                    end
                end
            end
  
            if tch
                for i = 1:length(Lambda)
                    for j = 1:n
                        if (s(j)==0 && W(j,Lambda(i))~=0)
                            a = W(j,Lambda(i))/t(Lambda(i));
                            if s(j)~=a
                                sc = true;
                                K(length(K)+1) = j;
                            end           
                            s(j) = a;
                        end
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (Wch || sch)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        %   recalculation of t 
            if Wch
                for i = 1:length(Omega)
                    i1 = Omega{i}(1);
                    i2 = Omega{i}(2);
                    if s(i1)~=0
                        a = W(i1,i2)/s(i1);
                        if t(i2)~=a
                            tc = true;
                            Lambda(length(Lambda)+1) = i2;
                        end 
                        t(i2) = a;
                    end
                end
            end
  
            if sch
                for i=1:length(K)
                    for j = 1:n
                        if (t(j)==0 && W(K(i),j)~=0)
                            a = W(K(i),j)/s(K(i));
                            if t(j)~=a
                                tc = true;
                                Lambda(length(Lambda)+1) = j;
                            end 
                            t(j) = a;
                        end
                    end
                end
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        end
        Wch = Wc;
        sch = sc;
        tch = tc;
        Wc = false;
        sc = false;
        tc = false;
  
    end 
end 
if gamma<omega 
    %   generation of new record of the list
    tL = length(L);  
    if (mode==1 && tL<szL) || (mode==2 && gamma<cstL)
        for l = 1:tL
            if L(l).gamma>gamma
                if l==1
                    l = 2;
                end
                L(l+1:tL+1) = L(l:tL);
                L = inlist(L,l,fl); 
                break
            end
            if l==tL 
                L = inlist(L,tL+1,fl); 
            end 
        end
    else
        if mode==1  
            for l = 1:szL
                if L(l).gamma>gamma
                    if l==1
                        l = 2;
                    end
                    c = length(Lk)+1;
                    Lk = transfer(Lk,L,c,szL);
                    L(l+1:szL) = L(l:szL-1);
                    L = inlist(L,l,fl);     
                    break 
                end
                if l==szL
                    Lk = inlist(Lk,length(Lk)+1,fl);
                end 
            end
        else
            Lk = inlist(Lk, length(Lk)+1,fl);
        end
    end
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function L = inlist(L,l,fl)  
%   adding new record in the list  
L(l).Q = Q;
L(l).r = r;
L(l).gamma = gamma;
L(l).x = vx;
L(l).W = W; 
L(l).s = s; 
L(l).t = t; 
if fl
    L(l).Y = inv_HBR(Q);
else
    L(l).Y = invQ;
end
end  
end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L1 = transfer(L1,L2,k1,k2)
     L1(k1).Q = L2(k2).Q;
     L1(k1).r = L2(k2).r;
     L1(k1).gamma = L2(k2).gamma;
     L1(k1).x = L2(k2).x;
     L1(k1).Y = L2(k2).Y;
     L1(k1).W = L2(k2).W;
     L1(k1).s = L2(k2).s;
     L1(k1).t = L2(k2).t;
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function vx = HBR(A,b) 
%   interval estimation based on Hansen-Bliek-Rohn's procedure 
n = dim(A);
C = inv(mid(A));
A = C*A;
b = C*b;
dA = diag(A); 
A = compmat(A);          %   comparison matrix  
B = inv(A);
v = abss(B*ones(n,1));
setround(-1)
u = A*v;
if ~all(min(u)>0) 
    error('Matrix of the system not an H-matrix')
else
    dAc = diag(A);
    A = A*B-eye(n); 
    setround(1)
    w = max(-A./(u*ones(1,n))); 
    dlow = v.*w'-diag(B);
    dlow = -dlow;
    B = B+v*w; 
    u = B*abss(b);
    d = diag(B);
    alpha = dAc+(-1)./d;
    k = size(b,2);
    if k==1
        beta = u./dlow - mag(b);
        vx = (b+midrad(0,beta))./(dA+midrad(0,alpha));
    else
        v = ones(1,k);
        beta = u./(d*v) - mag(b);
        vx = ( b + midrad(0,beta) ) ./ ( ( dA + midrad(0,alpha) ) * v );
    end
end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
function r = inv_HBR(a) 
%   computation of inverse interval matrix
  n = dim(a); 
  r = HBR(a,speye(n));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%