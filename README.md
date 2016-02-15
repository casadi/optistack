[![Build Status](https://travis-ci.org/casadi/optistack.png?branch=master)](https://travis-ci.org/casadi/optistack)

# optistack
The goal of this project is to provide a Yalmip-like wrapper around the Matlab interface of [CasADi](http://casadi.org),  resulting in an easy and fast way to solve (parametric) NLPs.

Example:
```matlab
x = optivar();
y = optivar();

nlp = optisolve((1-x)^2+100*(y-x^2)^2,{x^2+y^2<=1, x+y>=0});

optival(x)
optival(y)
```

Installation:
 * Obtain CasADi for your version of Matlab:

Windows   |   Linux     |    Mac
----------|-------------|--------------
[R2014a](http://files.casadi.org/3.0.0-rc2/windows/casadi-matlabR2014a-v3.0.0-rc2.zip) or later |    [R2014a](http://files.casadi.org/3.0.0-rc2/linux/casadi-matlabR2014a-v3.0.0-rc2.tar.gz) or later      | [R2015a](http://files.casadi.org/3.0.0-rc2/osx/casadi-matlabR2015a-v3.0.0-rc2.tar.gz) or later
[R2013a](http://files.casadi.org/3.0.0-rc2/windows/casadi-matlabR2013a-v3.0.0-rc2.zip) or later | | [R2014a](http://files.casadi.org/3.0.0-rc2/osx/casadi-matlabR2014a-v3.0.0-rc2.tar.gz) or later |

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
`optivar` inherits from `casadi.MX` class. This means [all the usual operations](http://casadi.sourceforge.net/v3.0.0-rc2/api/html/d9/dc2/group__expression__tools.html) for `MX` are possible.



Note: there is also a [python version](https://github.com/casadi/python-optistack/) available.

