function xf = rk4(ode,h,x,u)
  k1 = ode(x,       u);
  k2 = ode(x+h/2*k1,u);
  k3 = ode(x+h/2*k2,u);
  k4 = ode(x+h*k3,  u);
  xf = x + h/6 * (k1 + 2*k2 + 2*k3 + k4); 
end