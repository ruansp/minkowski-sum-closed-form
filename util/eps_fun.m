function val = eps_fun(x, epsilon)
%eps_fun exponentiation function
val = sign(x).*abs(x).^epsilon;
end