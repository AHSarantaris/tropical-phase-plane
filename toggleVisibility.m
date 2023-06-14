function toggleVisibility(varargin)
for i = 1:nargin
    varargin{i}.Visible = ~varargin{i}.Visible;
end
end

