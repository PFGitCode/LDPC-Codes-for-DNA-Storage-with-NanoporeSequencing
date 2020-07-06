function p = theta2p(fixPoint,theta,P)
if theta(fixPoint)~=0
    p= P*prod(theta(theta ~= 0))/theta(fixPoint);
else
    p= P*prod(theta(theta ~= 0));
end
% if abs(p) < 1.0e-5
%     p = 1.0e-10*(p/abs(p));
% end
end
