% 带约束最小化问题
function test
x0 = [0 0 0 0 0  0].';
a = [];
b = [];
lb = [0 0 0 0 0 0 ].';
ub = [inf inf inf inf inf inf].';
s = fmincon(@fun,x0,a,b,[],[],lb,ub,@con);
end


% 目标函数
function [y] = fun(X)
y =X(1)+X(2)+X(3)+X(4)+X(5)+X(6);
end

% 非线性约束
function [c,ceq]=cof(X)
A=100000;
c = A^6-(X(1)*A^5)+(A-X(1)*X(2)*A^4 + (A-X(1))*(A-X(2)) *A^3+(A- X(1))*(A-X(2))*(A-X(3))*X(4)*A^2 - (A-X(1))*(A-X(2))*(A-X(3))*(A-X(4))*X(5) *A+ (A-X(1))*(A-X(2))*(A-X(3))*(A-X(4))*(A-X(5))*X(6); % 造价-预算<=0
ceq = 0;
end


