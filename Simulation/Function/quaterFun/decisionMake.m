function result = decisionMake(theta,P)
    multiResult = P;
    for i = 1:length(theta)
         if theta(i) ~=0
             multiResult = multiResult*theta(i);
         end
    end
    result = multiResult;
end