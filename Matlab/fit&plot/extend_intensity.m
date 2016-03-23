% Extend and refine the x vector
%
%SYNOPSYS
% ex = EXTEND_INTENSITY(x)
%

function ex = extend_intensity(x)

n       = length(x);
m_x     = (max(x) - min(x))/2;
ex      = (linspace(max(0,min(x)-0.2*m_x), max(x)+0.2*m_x, 10*n)).';

end