# optistack
A simple matlab interface to casadi

The goal of this project is to provide a Yalmip-like Matlab interface to casadi.

```matlab
x = optivar()
y = optivar()
minimize((1-x)**2+100*(y-x**2)**2,[x**2+y**2<=1, x+y>=0])

value(x)
value(y)
```
