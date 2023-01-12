function nsqr=sellmeier(B,lambdar,lambda)
    nsqr = ones(size(lambda));
    for i=1:length(B)
        nsqr = nsqr + B(i).*lambda.^2./(lambda.^2-lambdar(i)^2);
    end
end