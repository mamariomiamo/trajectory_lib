%% evalute polynomial
% a is the derivative order, t is time to evaluate, order is polynomial order
function vec = poly_evaluate(a,t,order)
    vec = zeros(1,order+1);
    if a==0
        for i=1:order+1
            vec(i)=t^(i-1);
        end
    elseif a==1
        for i=2:order+1
            vec(i)= (i-1)*t^(i-2);
        end
    elseif a==2
        for i=3:order+1
            vec(i)= (i-2)*(i-1)*t^(i-3);
        end
    end
end