# optistack
The goal of this project is to provide a Yalmip-like wrapper around the Matlab interface of [casadi](http://casadi.org).

Example:
```matlab
x = optivar();
y = optivar();

nlp = optisolve((1-x)^2+100*(y-x^2)^2,{x^2+y^2<=1, x+y>=0});

optival(x)
optival(y)
```

Installation:
 * Obtain CasADi for [Windows](files.casadi.org/2.4.1/windows/casadi-matlabR2014a-v2.4.1.zip) or [Linux](files.casadi.org/2.4.1/linux/casadi-matlabR2014a-v2.4.1.tar.gz) and unzip it.
 * Add the unzipped directory to the Matlab path (`addpath('casadi_unzippeddir')`)
 * Obtain [Optistack](https://github.com/casadi/optistack/archive/master.zip) and unzip it.
 * Add the src directory to the Matlab path (`addpath('optistack_unzippeddir/src')`)
 * Run the above example to verify that the install succeeded


Cheat sheet:

                          |  Yalmip                     | Optistack
------------------------- | --------------------------- | -----------------------------
Declare decision variable | `x = sdpvar(nrows,ncols)`   | `x = optivar(nrows,ncols)`
Set initial values        | `assign(x,x0)`              | `x.setInit(x0)`
Solve a problem           | `optimize([x==0],f)`        | `optisolve(f,{x==0})`
Obtain numeric results    | `value(x)`                  | `optival(x)`


More functionality:
`optivar` inherits from `casadi.MX` class. This means [all the usual operations](http://casadi.sourceforge.net/v2.4.1/api/html/d9/dc2/group__expression__tools.html) for `MX` are possible.






