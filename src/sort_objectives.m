function [ scalar_objectives, twonorm_objectives, total_scalar_objective, total_objective ] = sort_objectives( objectives )

  if ~iscell(objectives)
    objectives = {objectives};
  end
  scalar_objectives = {};
  twonorm_objectives = {};
  total_objective = 0;
  total_scalar_objective = 0;
  for i=1:length(objectives)
    obj = objectives{i};
    if isvector(obj) && numel(obj)>1
      F = vec(obj);
      twonorm_objectives = {twonorm_objectives{:} F};
      total_objective = total_objective + 0.5*dot(F,F);
    else
      scalar_objectives = {scalar_objectives{:} obj};
      total_scalar_objective = total_scalar_objective + obj;
    end
  end
  total_objective = total_objective + total_scalar_objective; 
end

