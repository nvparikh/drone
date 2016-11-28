function Parameter = nondimensional_optimization(objectiveFun, Parameter_0, tMin, tMax, pos_tmin, pos_tmax, velocity_tmin, velocity_tmax, lb, ub)
   A_EQ = generate_A_EQ(tMin, tMax);
   B_eq = [pos_tmin;pos_tmax;velocity_tmin;velocity_tmax];
   integral_objectiveFun = @(tParameter) numericIntegration(@(t) objectiveFun(t, tParameter), tMin, tMax);
   Parameter = fmincon(integral_objectiveFun, Parameter_0, [], [], A_EQ, B_eq, lb, ub);
   t = tMin:0.05:tMax;
   position = fliplr(Parameter) * [t.^6; t.^5; t.^4; t.^3; t.^2; t; ones(1, size(t,2))];
   figure(1);
   plot(t, position);
   title('position');
   figure(2);
   velocity = fliplr(Parameter) * [6*t.^5; 5*t.^4; 4*t.^3; 3*t.^2; 2*t; ones(1, size(t,2)); zeros(1,size(t,2))];
   plot(t, velocity);
   title('velocity');
   figure(3);
   acceleration = fliplr(Parameter) * [5*6*t.^4; 4*5*t.^3; 3*4*t.^2; 2*3*t; 2*ones(1, size(t,2)); zeros(1,size(t,2)); zeros(1,size(t,2))];
   plot(t,acceleration);
   title('acceleration');
   figure(4)
   snap_square = objectiveFun(t, Parameter);
   plot(t, snap_square);
   title('snap_square');
end

function integration = numericIntegration(fun, tMin, tMax)
    integration = integral(fun, tMin, tMax);
end

function A_EQ = generate_A_EQ(tmin, tmax)
    A_EQ = [1 tmin tmin^2 tmin^3 tmin^4 tmin^5 tmin^6; 
            1 tmax tmax^2 tmax^3 tmax^4 tmax^5 tmax^6;
            0 1 2*tmin 3*tmin^2 4*tmin^3 5*tmin^4 6*tmin^5;
            0 1 2*tmax 3*tmax^2 4*tmax^3 5*tmax^4 6*tmax^5];
end

function val = test(t)
  val = [20*360 1*240 34*120 0 0 0 0] * [t.^2 t.^ t.^4 t.^3 t.^2 t 1]';
end
