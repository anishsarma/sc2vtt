function [R,M,J] =SynthSLS(T,A,B,J)
cvx_begin quiet
cvx_precision best
cvx_solver sedumi

if isnumeric(A)
    Atmp = A;
    A = {};
    for tt = 1:T
    A{tt} = Atmp;
    end
end
N = size(A{1},1);


if ~exist('J','var')
for tt = 1:T
    J.u{tt} = ones(N);
    J.x{tt} = ones(N);

end
J.uncertainEps = 0.0;
J.norm = inf;
end



if ~exist('B','var')
B = eye(N);
end
if isnumeric(B)
    Btmp = B;
    B = {};
    for tt = 1:T
        B{tt} = Btmp;
        if tt <= J.ifnDelay
            B{tt}(2,2) = 0;
        end
        if tt <= J.allDelay
            B{tt} = B{tt} * 0;
        end
    end
end
if J.ifnDelay <= 20
assignin('base','Bcell',B)
end
L = size(B{1},2);

variable R(N,N*T)
variable M(L,N*T)
variable t
Ju = zeros(L,N*T);
Jx = zeros(N,N*T);

cfx = 0;


t>=0
R(:,(1:N)) == t*eye(N);

R(:,(1:N)+(T-1)*N) == 0;
M(:,(1:N)+(T-1)*N) == 0;

R >= 0
M >= 0
for ii = 1:T
    for jj = 1:N
        if sum(B{ii}(jj,:)) == 0
            M(jj,(1:N)+(ii-1)*N) == 0;
        end
    end
end
for tt = 2:T
    R(:,(1:N)+(tt-1)*N) == A{tt-1} * R(:,(1:N)+(tt-2)*N) - B{tt-1} * M(:,(1:N)+(tt-2)*N);
end
for tt = 1:T
Ju(:,(1:N)+(tt-1)*N) = J.u{tt};
Jx(:,(1:N)+(tt-1)*N) = J.x{tt};
end

norm(J.uncertainEps*[R;M],J.norm)  <= t-1
t >= 1
cfx = norm([Jx(:,1:N:end).*R(:,1:N:end);Ju(:,1:N:end).*M(:,1:N:end)],J.norm);
cfx = norm([Jx.*R;Ju.*M],J.norm);

minimize(cfx)

J.cfx = cfx;
J.umat = Ju;
J.xmat = Jx; 
cvx_end

R = full(R/t);
M = full(M/t);