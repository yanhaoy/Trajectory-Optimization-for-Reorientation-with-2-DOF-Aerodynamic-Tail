function [out] = quat2rotm(in)

[qr, qi, qj, qk] = deal(in(1), in(2), in(3), in(4));

out = [1-2*(qj^2+qk^2), 2*(qi*qj-qk*qr), 2*(qi*qk+qj*qr);
    2*(qi*qj+qk*qr), 1-2*(qi^2+qk^2), 2*(qj*qk-qi*qr);
    2*(qi*qk-qj*qr), 2*(qj*qk+qi*qr), 1-2*(qi^2+qj^2)];

end

