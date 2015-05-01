function [ out ] = poor_veccat( varargin )

    out_cell = {};
    for i=1:length(varargin)
        out_cell = {out_cell{:},vec(varargin{i})};
        
    end
    out = vertcat(out_cell{:});
end

