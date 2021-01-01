myfun = @(z,omega) z*sin(2*omega) + sin(2*omega*z);  % parameterized function
omega = pi*3/4;
fun = @(z) myfun(z,omega);    % function of z alone
format long
z = fzero(fun,0.5)
fun(z)
format