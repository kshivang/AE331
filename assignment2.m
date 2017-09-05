function assignment2()
 

%Part A
disp('Part a:')
for i = 1:5
    for j = 1:5
        for k = 1:5
            [V, D] = eig(getRho(0.2 * i, 0.1 * j, 0.05 * k));
            print = ['For i = ', num2str(i), ', j = ', num2str(j), ', k = ', num2str(k)];
            disp(print)
            disp('Eigen Value matrix is')
            disp(D)
            disp('Eigen Vector matrix is')
            disp(V)
        end
    end
end



% Part B
disp('Part b:')
[V, D] = eig(getRho(1, 1, 1));
disp('For Point(1, 1, 1)')
disp('Eigen Value matrix is')
disp(D)
disp('Eigen Vector(V) matrix is')
disp(V)

checkOrthogonality(V);
disp('Orthogonal Eigen Vectors')


OSS = getOCC(D(1, 1), D(2, 2), D(3, 3));
disp('Octahedral Shear Shtress at 1, 1, 1 is')
disp(OSS)

ShearMax = getShearMax(D(1, 1), D(2, 2), D(3, 3));
disp('Shear Max at 1, 1, 1 is')
disp(ShearMax)


[V, D] = eig(getRho(5, 5, 5));
disp('For Point(5, 5, 5)')
disp('Eigen Value matrix is')
disp(D)
disp('Eigen Vector matrix is')
disp(V)

checkOrthogonality(V);
disp('Orthogonal Eigen Vectors')

OSS = getOCC(D(1, 1), D(2, 2), D(3, 3));
disp('Octahedral Shear Shtress at 5, 5, 5 is')
disp(OSS)

ShearMax = getShearMax(D(1, 1), D(2, 2), D(3, 3));
disp('Shear Max at 5, 5, 5 is')
disp(ShearMax)



% Part C
rho = getRho(3, 3,3);

n = [(1/3)^(1/2),  (1/3)^(1/2), (1/3)^(1/2)];
StressVector = (rho * (n.')).';
disp('Stress Vector for n')
disp(StressVector)

p = dot(StressVector, n) * n;
disp('Stress Vector component perpendicular to Plane n')
disp(p)

pn = StressVector - p;
disp('Stress Vector component parallel to Plane(n)')
disp(pn)


m = [(3^(1/2))/2,  (1/2), 0];
StressVector = (rho * (n.')).';
disp('Stress Vector for m')
disp(StressVector)

p = dot(StressVector, m) * m;
disp('Stress Vector component perpendicular to Plane m')
disp(p)

pm = StressVector - p;
disp('Stress Vector component parallel to Plane(m)')
disp(pm)


% Part D
rho = getRho(3, 3, 3);
rotatedRho = rotatateTensorInX(rho, 30);
disp('Stress Vector in rotated coordinate system is given by')
disp(rotatedRho)

% Part E
rho = getRho(3, 3, 3);
[V, D] = eig(rho);
disp('For Point(3, 3, 3)')
disp('Eigen Value matrix is')
disp(D)
disp('Eigen Vector(V) matrix is')
disp(V)
rotatedRho = rotatateTensorInX(rho, 30);
disp('Stress Vector in rotated coordinate system is given by')
disp(rotatedRho)
[V, D] = eig(rho);
disp('For Point(3, 3, 3) rotated by 30degree about X')
disp('Eigen Value matrix is')
disp(D)
disp('Eigen Vector(V) matrix for rotated coodinate is')
disp(V)
disp('Therefor eigenvalues are actually independent of coordinate rotation')

    %%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%
    
    function rho = getRho(x1, x2, x3)
        rho = ...
            [ ((10*x1^2) - (3 * x2))                (-3*x1*x2)               (8*x1*x3); ...
              (-3*x1*x2)                                  (-5*x2*x3)               (3*x2^2);  ...
              (8*x1*x3)                                   (3*x2^2)                 (-10*x3)];
    end

    function OCC = getOCC(lemda1, lemda2, lemda3)
        OCC = (((lemda1 - lemda2)^2 + (lemda2 - lemda3)^2 + (lemda3 - lemda1)^2)^(1/2)) / 3;
    end

    function ShearMax = getShearMax(rhoX, rhoY, rhoZ)
        SMax(1) = abs((rhoX - rhoY) / 2);
        SMax(3) = abs((rhoY - rhoZ) / 2);
        SMax(3) = abs((rhoX - rhoZ) / 2);
        ShearMax = max(SMax);
    end
    
    function isOrthogonal = checkOrthogonality(V)
        VTranspose = V.';
        disp('Transpose of V ')
        disp(VTranspose)
        
        OrthogonalityCheck = V*VTranspose;
        disp('V * VTranspose')
        disp(OrthogonalityCheck)
        isOrthogonal = ((OrthogonalityCheck - eye(3)) == zeros(3));
    end

    function rotatedTensorX = rotatateTensorInX(tensor, degree)
        RotationMatrix = [1                             0                         0 ;...
                                     0                    cosd(degree)          -sind(degree);...
                                     0                    sind(degree)            cosd(degree)];
         rotatedTensorX = (RotationMatrix.') * tensor * RotationMatrix;
    end
end
