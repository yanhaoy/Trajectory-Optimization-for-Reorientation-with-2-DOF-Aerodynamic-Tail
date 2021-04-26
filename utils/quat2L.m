function [out] = quat2L(in)

[e0, e1, e2, e3] = deal(in(1), in(2), in(3), in(4));

out = [-e1, e0, e3, -e2;
    -e2, -e3, e0, e1;
    -e3, e2, -e1, e0];

end

