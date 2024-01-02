% MATLAB Codes for 4Node Quadrilateral Element
% © Hamed Valaei


% Clear workspace
clear all
clc

% Load data from text files
Nodes = load('Nodes.txt');
Forces = load('Forces.txt');
Elements = load('Elements.txt');

% Constrained Node and Direction
[ConstrainedNodes, ConstrainedNodeDirection] = textread('Constraints.txt', '%d %c');

% Sizes
[NumNodes, ~] = size(Nodes);
[NumElements, ~] = size(Elements);
[NumForces, ~] = size(Forces);
[NumConstraints, ~] = size(ConstrainedNodes);

% Allocate coordinates of nodes to X and Y vectors
X = Nodes(:, 2);
Y = Nodes(:, 3);

% Material properties
E = input('Enter Modulus of Elasticity (E): ');
p = input('Use p = 1 for cases of plane stress, and p = 2 for cases of plane strain: ');
nu = input('Input Poisson’s ratio: ');
t = input('Input Thickness: ');

% Force vector
ForceVector = zeros(2 * NumNodes, 1);

% Populate force vector
for q = 1:NumForces
    % Locating position of xForces
    ForceVector(2 * Forces(q, 1) - 1) = Forces(q, 2);
    % Locating position of yForces
    ForceVector(2 * Forces(q, 1)) = Forces(q, 3);
end

% Allocate nodes of the element to i, j, m, and n
i = Elements(:, 2);
j = Elements(:, 3);
m = Elements(:, 4);
n = Elements(:, 5);

% Calculate D according to plane stress or plane strain
if p == 1
    D = (E / (1 - nu^2)) * [1 nu 0; nu 1 0; 0 0 (1 - nu) / 2];
elseif p == 2
    D = (E / (1 + nu) / (1 - 2 * nu)) * [1 - nu nu 0; nu 1 - nu 0; 0 0 (1 - 2 * nu) / 2];
end

syms r s;

% Shape Functions of Nodes
N(1) = 1/4 * (1 - r) * (1 - s);
N(2) = 1/4 * (1 + r) * (1 - s);
N(3) = 1/4 * (1 + r) * (1 + s);
N(4) = 1/4 * (1 - r) * (1 + s);

% Initialize matrices
KG = zeros(2 * NumNodes, 2 * NumNodes);
k = zeros(8, 8);

% Loop over elements
for q = 1:NumElements
    % Compute element properties
    a = 1/4 * (Y(i(q)) * (r - 1) + Y(j(q)) * (-1 - r) + Y(m(q)) * (1 + r) + Y(n(q)) * (1 - r));
    b = 1/4 * (Y(i(q)) * (s - 1) + Y(j(q)) * (1 - s) + Y(m(q)) * (1 + s) + Y(n(q)) * (-1 - s));
    c = 1/4 * (X(i(q)) * (s - 1) + X(j(q)) * (1 - s) + X(m(q)) * (1 + s) + X(n(q)) * (-1 - s));
    d = 1/4 * (X(i(q)) * (r - 1) + X(j(q)) * (-1 - r) + X(m(q)) * (1 + r) + X(n(q)) * (1 - r));

    % Determine B1, B2, B3, B4
    for qq = 1:4
        ds = diff(N(qq), s);
        dr = diff(N(qq), r);
        Bi(:, :, qq) = [a * dr - b * ds, 0; 0, c * ds - d * dr; c * ds - d * dr, a * dr - b * ds];
    end

    % Calculate Determinant of Jacobian
    deteJ = (1/8) * [X(i(q)), X(j(q)), X(m(q)), X(n(q))] * [0, 1 - s, s - r, r - 1; s - 1, 0, r + 1, -r - s;
        r - s, -r - 1, 0, s + 1; 1 - r, r + s, -s - 1, 0] * [Y(i(q)); Y(j(q)); Y(m(q)); Y(n(q))];

    % Final B Matrix
    BB = 1/deteJ * [Bi(:, :, 1), Bi(:, :, 2), Bi(:, :, 3), Bi(:, :, 4)];
    B(:, :, q) = BB;
    B2 = deteJ * B(:, :, q)' * D * B(:, :, q);

    % MATLAB Calculate integral for stiffness matrix
    integ = int(int(B2, s, -1, 1), r, -1, 1);
    k(:, :, q) = double(t * integ);

    % Assemble stiffness matrix of elements in Global Stiffness matrix
    KG(2 * i(q) - 1:2 * i(q), 2 * i(q) - 1:2 * i(q)) = KG(2 * i(q) - 1:2 * i(q), 2 * i(q) - 1:2 * i(q)) + k(1:2, 1:2, q);
    KG(2 * i(q) - 1:2 * i(q), 2 * j(q) - 1:2 * j(q)) = KG(2 * i(q) - 1:2 * i(q), 2 * j(q) - 1:2 * j(q)) + k(1:2, 3:4, q);
    KG(2 * i(q) - 1:2 * i(q), 2 * m(q) - 1:2 * m(q)) = KG(2 * i(q) - 1:2 * i(q), 2 * m(q) - 1:2 * m(q)) + k(1:2, 5:6, q);
    KG(2 * i(q) - 1:2 * i(q), 2 * n(q) - 1:2 * n(q)) = KG(2 * i(q) - 1:2 * i(q), 2 * n(q) - 1:2 * n(q)) + k(1:2, 7:8, q);
    KG(2 * j(q) - 1:2 * j(q), 2 * i(q) - 1:2 * i(q)) = KG(2 * j(q) - 1:2 * j(q), 2 * i(q) - 1:2 * i(q)) + k(3:4, 1:2, q);
    KG(2 * j(q) - 1:2 * j(q), 2 * j(q) - 1:2 * j(q)) = KG(2 * j(q) - 1:2 * j(q), 2 * j(q) - 1:2 * j(q)) + k(3:4, 3:4, q);
    KG(2 * j(q) - 1:2 * j(q), 2 * m(q) - 1:2 * m(q)) = KG(2 * j(q) - 1:2 * j(q), 2 * m(q) - 1:2 * m(q)) + k(3:4, 5:6, q);
    KG(2 * j(q) - 1:2 * j(q), 2 * n(q) - 1:2 * n(q)) = KG(2 * j(q) - 1:2 * j(q), 2 * n(q) - 1:2 * n(q)) + k(3:4, 7:8, q);
    KG(2 * m(q) - 1:2 * m(q), 2 * i(q) - 1:2 * i(q)) = KG(2 * m(q) - 1:2 * m(q), 2 * i(q) - 1:2 * i(q)) + k(5:6, 1:2, q);
    KG(2 * m(q) - 1:2 * m(q), 2 * j(q) - 1:2 * j(q)) = KG(2 * m(q) - 1:2 * m(q), 2 * j(q) - 1:2 * j(q)) + k(5:6, 3:4, q);
    KG(2 * m(q) - 1:2 * m(q), 2 * m(q) - 1:2 * m(q)) = KG(2 * m(q) - 1:2 * m(q), 2 * m(q) - 1:2 * m(q)) + k(5:6, 5:6, q);
    KG(2 * m(q) - 1:2 * m(q), 2 * n(q) - 1:2 * n(q)) = KG(2 * m(q) - 1:2 * m(q), 2 * n(q) - 1:2 * n(q)) + k(5:6, 7:8, q);
    KG(2 * n(q) - 1:2 * n(q), 2 * i(q) - 1:2 * i(q)) = KG(2 * n(q) - 1:2 * n(q), 2 * i(q) - 1:2 * i(q)) + k(7:8, 1:2, q);
    KG(2 * n(q) - 1:2 * n(q), 2 * j(q) - 1:2 * j(q)) = KG(2 * n(q) - 1:2 * n(q), 2 * j(q) - 1:2 * j(q)) + k(7:8, 3:4, q);
    KG(2 * n(q) - 1:2 * n(q), 2 * m(q) - 1:2 * m(q)) = KG(2 * n(q) - 1:2 * n(q), 2 * m(q) - 1:2 * m(q)) + k(7:8, 5:6, q);
    KG(2 * n(q) - 1:2 * n(q), 2 * n(q) - 1:2 * n(q)) = KG(2 * n(q) - 1:2 * n(q), 2 * n(q) - 1:2 * n(q)) + k(7:8, 7:8, q);
end

% Save a copy of the original Global Stiffness Matrix
KG2 = KG;

% Eliminate rows and columns of KG matrix according to boundary conditions
for r = 1:NC
    nn = 2 * CN(r, 1);
    % Nodes constrained in X Direction
    if CND(r, 1) == 'X'
        KG2(nn - 1, :) = 0;
        KG2(:, nn - 1) = 0;
        KG2(nn - 1, nn - 1) = 1;
    % Nodes constrained in Y Direction
    elseif CND(r, 1) == 'Y'
        KG2(nn, :) = 0;
        KG2(:, nn) = 0;
        KG2(nn, nn) = 1;
    % Nodes constrained in Both X and Y Directions
    elseif CND(r, 1) == 'B'
        KG2(nn - 1, :) = 0;
        KG2(:, nn - 1) = 0;
        KG2(nn - 1, nn - 1) = 1;
        KG2(nn, :) = 0;
        KG2(:, nn) = 0;
        KG2(nn, nn) = 1;
    end
end

% Compute displacement matrix
d = KG2 \ F;

% Apply constrained values to displacement matrix
for r = 1:NC
    nn = 2 * CN(r, 1);
    % Nodes constrained in X Direction
    if CND(r, 1) == 'X'
        d(nn - 1, 1) = 0;
    % Nodes constrained in Y Direction
    elseif CND(r, 1) == 'Y'
        d(nn, 1) = 0;
    % Nodes constrained in Both X and Y Directions
    elseif CND(r, 1) == 'B'
        d(nn - 1, 1) = 0;
        d(nn, 1) = 0;
    end
end

% Compute Support Reactions
R = KG * d;

% Initialize matrices for stress, strain, and principal stresses
sigma = zeros(3, NE);
epsilon = zeros(3, NE);
s1 = zeros(NE, 1);
s2 = zeros(NE, 1);
theta = zeros(NE, 1);

% Loop over elements to calculate stress, strain, and principal stresses
u = zeros(8, 1);
for q = 1:NE
    % Calculate u matrix for each element
    u = [d(2 * i(q) - 1); d(2 * i(q)); d(2 * j(q) - 1); d(2 * j(q)); d(2 * m(q) - 1); d(2 * m(q)); d(2 * n(q) - 1); d(2 * n(q))];
    % Calculate Stress at the centroid of elements
    sig = D * B(:, :, q) * u;
    sig1 = subs(sig, s, 0);
    sig2 = subs(sig1, r, 0);
    sigma(:, q) = double(sig2);
    epsilon(:, q) = sigma(:, q) \ D;
    U = (sigma(1, q) + sigma(2, q)) / 2;
    Q = ((sigma(1, q) - sigma(2, q)) / 2)^2 + sigma(3, q) * sigma(3, q);
    M = 2 * sigma(3, q) / (sigma(1, q) - sigma(2, q));
    % Principal stresses
    s1(q) = U + sqrt(Q);
    s2(q) = U - sqrt(Q);
    theta(q) = (atan(M) / 2) * 180 / pi;
end

% Display analysis results
disp('|/////////////////////////////////////////////////////////////////|')
disp('|/////// Analysis Results of 4Node Quadrilateral Element /////////|')
disp('|////////////© Hamed Valaei////////////////////////////////////////|')
disp('>>>>>>>> Nodal Displacement')
disp('__Node_________________dx_______________________dy__________________')
disp('   ')
for i = 1:NN
    fprintf(' %g\t', i);
    fprintf('                  %g\t', d(2 * i - 1));
    fprintf('       %g\t', d(2 * i));
    disp(' ')
end

disp('   ')
disp('>>>>>>>> Support Reactions')
disp('   ')
disp('__Node_________________Rx_______________________Ry__________________')
for q = 1:NC
    if CND(q, 1) == 'X'
        fprintf(' %g\t', CN(q, 1));
        fprintf('                  %g\t', R(2 * CN(q, 1) - 1))
        disp(' ')
    elseif CND(q, 1) == 'Y'
        fprintf(' %g\t', CN(q, 1));
        fprintf('                                                         %g\t', R(2 * CN(q, 1)))
        disp(' ')
    elseif CND(q, 1) == 'B'
        fprintf(' %g\t', CN(q, 1));
        fprintf('                  %g\t', R(2 * CN(q, 1) - 1));
        fprintf('                   %g\t', R(2 * CN(q, 1)))
        disp(' ')
    end
end

disp('   ')
disp('>>>>>>>> Linear Triangle Element Stresses')
disp('   ')
disp('__Element No.________Sigmax___________Sigmay__________Txy__________')
for q = 1:NE
    fprintf('   %g\t', q);
    fprintf('     %g\t', sigma(1, q));
    fprintf('      %g\t', sigma(2, q));
    fprintf('      %g\t', sigma(3, q));
    disp(' ')
end

disp('   ')
disp('>>>>>>>> Linear Triangle Element Strains')
disp('   ')
disp('__Element No.________epsx___________epsy__________gammaxy__________')
for q = 1:NE
    fprintf('   %g\t', q);
    fprintf('     %g\t', epsilon(1, q));
    fprintf('      %g\t', epsilon(2, q));
    fprintf('      %g\t', epsilon(3, q));
    disp(' ')
end

disp('   ')
disp('>>>>>>>> Linear Triangle Element Principal Stresses')
disp('   ')
disp('__Element No.________Sig1___________Sig2__________Tmax__________')
for q = 1:NE
    fprintf('   %g\t', q);
    fprintf('     %g\t', s1(q));
    fprintf('      %g\t', s2(q));
    fprintf('      %g\t', theta(q));
    disp(' ')
end


