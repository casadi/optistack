function [] = callback(fun)
  global x
  global y
  plot(optival(x),optival(y),'ro')
  pause(1)
end
