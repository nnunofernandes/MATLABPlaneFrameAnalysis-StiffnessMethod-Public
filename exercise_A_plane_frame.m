%% EXERCISE A -- PLANE FRAME BY THE DIRECT STIFFNESS METHOD
% Standalone teaching script. Run the complete file, or run each %% section.
% Local functions are included below; MATLAB R2016b or later is required.
%
% MODEL: linear elastic, small-displacement Euler-Bernoulli frame elements.
% Each node has three DOFs: ux, uy, theta. Axial deformation is included;
% shear deformation and geometric nonlinearity are excluded.
%
% UNITS: N, mm, rad. Moments are N*mm internally.
% GLOBAL SIGNS: +x right, +y up, +theta counterclockwise (CCW).
% The DOF NUMBERS follow the supplied figure. The rotation SIGN convention
% is explicitly CCW throughout this script, regardless of drawn arrows.
% LOCAL AXES: x from element node i to j; y is 90 degrees CCW from local x.
% Local vector order: [ui vi thetai uj vj thetaj].
%
% Main workflow:
% geometry -> element stiffness -> coordinate transformation -> assembly
% -> boundary conditions -> solve -> reactions/end forces -> verification.

clear; clc; close all;
format short g;
showMatrices = true;   % Show element matrices and assembled equations.
makePlots = true;      % Set false for numerical calculations only.
deformationScale = 0;  % Zero selects an automatic visual magnification.

%% 1. Input data, geometry, connectivity and DOF numbering
% Node 1: loaded free end; node 2: rigid corner; node 3: fixed base.
% Element 1 runs 1 -> 2 (right); element 2 runs 2 -> 3 (down).
E = 200e3;                 % N/mm^2 = 200 GPa
I = 300e6;                 % mm^4
A = 12.9e3;                % mm^2
nodes = [0 4000;           % [x y], mm
         3000 4000;
         3000 0];
elements = [1 2; 2 3];
dofMap = [1 2 3;           % Node 1: ux uy theta
          4 5 6;           % Node 2
          7 8 9];          % Node 3
fixedDofs = [7 8 9];       % Fixed support: both translations and rotation.
prescribedValues = [0; 0; 0];

% POINT LOADS: 18 kN left and 27 kN down at node 1; no applied couple.
Fpoint = zeros(9,1);
Fpoint(dofMap(1,1)) = -18e3;
Fpoint(dofMap(1,2)) = -27e3;

% Uniform distributed loads in ELEMENT LOCAL axes: [qx qy], N/mm.
qLocal = zeros(2,2);       % No distributed loads in this exercise.
fprintf('EXERCISE A: displacements at node 1\n');
fprintf('Global signs: right/up/CCW positive. Internal units: N, mm, rad.\n');

%% 2. Element stiffness, coordinate transformation and load assembly
nNode = size(nodes,1);
nElem = size(elements,1);
nDof = 3*nNode;
assert(isequal(sort(dofMap(:)), (1:nDof).'), 'Check DOF numbering.');
K = zeros(nDof);
F = Fpoint;               % Total equivalent GLOBAL nodal load vector.
elementData = struct([]);

for e = 1:nElem
    ij = elements(e,:);
    delta = nodes(ij(2),:) - nodes(ij(1),:);
    L = hypot(delta(1),delta(2));
    c = delta(1)/L;
    s = delta(2)/L;
    edofs = [dofMap(ij(1),:) dofMap(ij(2),:)];

    % q_global = T' * q_local and d_local = T * d_global.
    % Therefore k_global = T' * k_local * T.
    [kLocal,T] = frameElement(E,A,I,L,c,s);
    kGlobal = T.'*kLocal*T;

    % CONSISTENT load vector for constant local qx and qy:
    % [qx*L/2; qy*L/2; qy*L^2/12; qx*L/2; qy*L/2; -qy*L^2/12].
    % The moments are essential: two equal transverse forces alone are not
    % the work-equivalent loading of a bending element.
    qx = qLocal(e,1);
    qy = qLocal(e,2);
    fLocal = [qx*L/2; qy*L/2; qy*L^2/12;
              qx*L/2; qy*L/2; -qy*L^2/12];
    fGlobal = T.'*fLocal;

    K(edofs,edofs) = K(edofs,edofs) + kGlobal;
    F(edofs) = F(edofs) + fGlobal;

    % Retain each element's data for force recovery and plotting.
    elementData(e).L = L;
    elementData(e).c = c;
    elementData(e).s = s;
    elementData(e).dofs = edofs;
    elementData(e).kLocal = kLocal;
    elementData(e).T = T;
    elementData(e).fLocal = fLocal;

    fprintf('\nElement %d: node %d -> %d, L = %.0f mm, c = %.1f, s = %.1f\n', ...
        e,ij(1),ij(2),L,c,s);
    fprintf('Element DOFs: '); fprintf('%d ',edofs); fprintf('\n');
    if showMatrices
        disp('Local stiffness kLocal:'); disp(kLocal);
        disp('Transformation T (dLocal = T*dGlobal):'); disp(T);
        disp('Global element stiffness T''*kLocal*T:'); disp(kGlobal);
        disp('Equivalent local nodal load [N N N*mm N N N*mm]'':');
        disp(fLocal);
    end
end

%% 3. Apply support conditions and solve the reduced system
% Partition: Kff*df + Kfc*dc = Ff, where dc contains known displacements.
% Backslash solves the system directly; an explicit inverse is unnecessary.
freeDofs = setdiff(1:nDof,fixedDofs);
d = zeros(nDof,1);
d(fixedDofs) = prescribedValues;
Kff = K(freeDofs,freeDofs);
rhs = F(freeDofs) - K(freeDofs,fixedDofs)*d(fixedDofs);

assert(norm(K-K.','fro') <= 1e-12*norm(K,'fro'), 'K must be symmetric.');
[~,cholFlag] = chol(Kff);
assert(cholFlag == 0, 'The restrained frame is unstable or incorrectly assembled.');
d(freeDofs) = Kff \ rhs;

if showMatrices
    fprintf('\nASSEMBLED SYSTEM\n');
    disp('Full K (mixed force/moment and translation/rotation units):'); disp(K);
    disp('Equivalent global load vector F:'); disp(F);
    disp('Free DOFs, in reduced-system order:'); disp(freeDofs);
    disp('Reduced stiffness Kff:'); disp(Kff);
    disp('Reduced right-hand side:'); disp(rhs);
end

%% 4. Recover reactions and member-end actions
% R = K*d-F includes distributed loads through their assembled equivalent
% nodal vector. Free-DOF entries should be numerical roundoff only.
R = K*d - F;
nodalD = reshape(d(dofMap.'),3,nNode).';
nodalR = reshape(R(dofMap.'),3,nNode).';
displacementTable = table((1:nNode).',dofMap(:,1),dofMap(:,2),dofMap(:,3), ...
    nodalD(:,1),nodalD(:,2),nodalD(:,3), ...
    'VariableNames',{'Node','DOF_ux','DOF_uy','DOF_theta','ux_mm','uy_mm','theta_rad'});
fprintf('\nNODAL DISPLACEMENTS (positive right/up/CCW)\n');
disp(displacementTable);

supportNodes = find(any(ismember(dofMap,fixedDofs),2));
reactionTable = table(supportNodes,nodalR(supportNodes,1)/1e3, ...
    nodalR(supportNodes,2)/1e3,nodalR(supportNodes,3)/1e6, ...
    'VariableNames',{'Node','Rx_kN','Ry_kN','Mz_kNm'});
fprintf('SUPPORT REACTIONS (positive right/up/CCW)\n');
disp(reactionTable);

endForcesLocal = zeros(nElem,6);
for e = 1:nElem
    data = elementData(e);
    dLocal = data.T*d(data.dofs);
    % Subtract the consistent load vector. kLocal*dLocal alone gives the
    % wrong end actions when the member carries distributed loading.
    qEnd = data.kLocal*dLocal - data.fLocal;
    elementData(e).dLocal = dLocal;
    endForcesLocal(e,:) = qEnd.';
end
endForceTable = table((1:nElem).',endForcesLocal(:,1)/1e3, ...
    endForcesLocal(:,2)/1e3,endForcesLocal(:,3)/1e6, ...
    endForcesLocal(:,4)/1e3,endForcesLocal(:,5)/1e3,endForcesLocal(:,6)/1e6, ...
    'VariableNames',{'Element','Ni_kN','Vi_kN','Mi_kNm','Nj_kN','Vj_kN','Mj_kNm'});
fprintf('LOCAL END ACTIONS applied TO each member, in its local axes\n');
fprintf('Order: [Ni Vi Mi Nj Vj Mj]; moments positive CCW at both ends.\n');
fprintf('These are end actions, not a bending-diagram sign convention.\n');
disp(endForceTable);

%% 5. Numerical and overall static-equilibrium checks
% Sum the physical forces/moments (distributed loads at their centroids).
% This independently checks the assembled consistent-load vector.
pointNodal = reshape(Fpoint(dofMap.'),3,nNode).';
balance = [sum(pointNodal(:,1)+nodalR(:,1));
           sum(pointNodal(:,2)+nodalR(:,2));
           sum(pointNodal(:,3)+nodalR(:,3) ...
             + nodes(:,1).*(pointNodal(:,2)+nodalR(:,2)) ...
             - nodes(:,2).*(pointNodal(:,1)+nodalR(:,1)))];
for e = 1:nElem
    data = elementData(e);
    xy = mean(nodes(elements(e,:),:),1);
    qGlobal = [data.c -data.s; data.s data.c]*qLocal(e,:).';
    resultant = qGlobal*data.L;
    balance = balance + [resultant;
        xy(1)*resultant(2)-xy(2)*resultant(1)];
end
forceScale = max(1,norm(F,inf));
lengthScale = max([elementData.L]);
assert(norm(R(freeDofs),inf) < 1e-9*forceScale, 'Free-DOF equilibrium failed.');
assert(norm(balance(1:2),inf) < 1e-8*forceScale, 'Force equilibrium failed.');
assert(abs(balance(3)) < 1e-8*forceScale*lengthScale, 'Moment equilibrium failed.');
fprintf('Equilibrium residual: Fx = %.3e N, Fy = %.3e N, M = %.3e N*mm\n',balance);

%% 6. Independent analytical check: virtual work / cantilever formulas
% Let P and Q be the positive magnitudes of the leftward/downward loads.
% The column-top moment caused by Q is Q*B, and the top rotation is
% theta2 = (P*H^2/2 + Q*B*H)/(E*I), positive CCW.
% ux1 = column bending displacement + beam axial shortening.
% uy1 = column shortening + corner-rotation contribution + beam bending.
B = 3000; H = 4000; P = 18e3; Q = 27e3;
theta2 = (P*H^2/2 + Q*B*H)/(E*I);
d1Analytical = [-P*H^3/(3*E*I) - Q*B*H^2/(2*E*I) - P*B/(E*A);
               -Q*H/(E*A) - B*theta2 - Q*B^3/(3*E*I);
               theta2 + Q*B^2/(2*E*I)];
d1 = d(dofMap(1,:));
assert(all(abs(d1-d1Analytical) < 1e-8*max(abs(d1Analytical),1e-6)), ...
    'Independent displacement check failed.');
R3Analytical = [P; Q; -(P*H+Q*B)];
assert(all(abs(R(dofMap(3,:))-R3Analytical) < 1e-7*max(abs(R3Analytical),1)), ...
    'Independent support-reaction check failed.');
fprintf('\nREQUESTED ANSWER -- NODE 1\n');
fprintf('ux1    = %+.8f mm  (negative = left)\n',d1(1));
fprintf('uy1    = %+.8f mm  (negative = down)\n',d1(2));
fprintf('theta1 = %+.9f rad = %+.6f deg (positive = CCW)\n',d1(3),d1(3)*180/pi);
fprintf('Analytical and stiffness-method results agree.\n');

%% 7. Plot undeformed and deformed shapes
if makePlots
    plotFrame(nodes,elements,dofMap,d,elementData,qLocal,E,A,I, ...
        supportNodes,fixedDofs,deformationScale,'Exercise A');
    subplot(1,2,1); hold on;
    quiver(0,4,-0.55,0,0,'r','LineWidth',1.8,'MaxHeadSize',0.45);
    quiver(0,4,0,-0.65,0,'r','LineWidth',1.8,'MaxHeadSize',0.45);
    text(-0.62,4.15,'18 kN','Color','r');
    text(0.1,3.55,'27 kN','Color','r');
end
% Useful workspace variables: K, F, Kff, d, R, elementData,
% displacementTable, reactionTable, endForceTable, d1Analytical.

%% LOCAL FUNCTIONS -- included so this script is fully standalone
function [k,T] = frameElement(E,A,I,L,c,s)
% 2-node, 6-DOF Euler-Bernoulli plane-frame stiffness in LOCAL axes.
% Axial stiffness: EA/L. Bending stiffness: EI times powers of 1/L.
    a = E*A/L;
    b = 12*E*I/L^3;
    g = 6*E*I/L^2;
    h = 4*E*I/L;
    j = 2*E*I/L;
    k = [ a  0  0 -a  0  0;
          0  b  g  0 -b  g;
          0  g  h  0 -g  j;
         -a  0  0  a  0  0;
          0 -b -g  0  b -g;
          0  g  j  0 -g  h];
    T = [ c  s  0  0  0  0;
         -s  c  0  0  0  0;
          0  0  1  0  0  0;
          0  0  0  c  s  0;
          0  0  0 -s  c  0;
          0  0  0  0  0  1];
end

function plotFrame(nodes,elements,dofMap,d,data,qLocal,E,A,I, ...
    supportNodes,fixedDofs,scale,caption)
% Bending uses Hermite interpolation plus the exact UDL particular solution.
% A straight line between displaced nodes would hide member curvature.
% The particular terms vanish at the ends (and their bending slopes vanish).
    nElem = size(elements,1);
    shape = cell(nElem,1);
    maxDisplacement = 0;
    for e = 1:nElem
        L = data(e).L;
        x = linspace(0,L,81);
        z = x/L;
        dl = data(e).dLocal;
        u = (1-z)*dl(1)+z*dl(4) + qLocal(e,1)*x.*(L-x)/(2*E*A);
        N1 = 1-3*z.^2+2*z.^3;
        N2 = L*(z-2*z.^2+z.^3);
        N3 = 3*z.^2-2*z.^3;
        N4 = L*(-z.^2+z.^3);
        v = N1*dl(2)+N2*dl(3)+N3*dl(5)+N4*dl(6) ...
            + qLocal(e,2)*x.^2.*(L-x).^2/(24*E*I);
        displacement = [data(e).c -data(e).s; data(e).s data(e).c]*[u;v];
        origin = nodes(elements(e,1),:).';
        xy = repmat(origin,1,numel(x)) + [data(e).c;data(e).s]*x;
        shape{e} = struct('xy',xy,'displacement',displacement);
        maxDisplacement = max(maxDisplacement,max(sqrt(sum(displacement.^2,1))));
    end
    span = max(max(nodes,[],1)-min(nodes,[],1));
    if scale <= 0
        scale = 0.12*span/max(maxDisplacement,eps);
    end
    figure('Color','w','Name',caption,'Position',[100 100 1100 520]);
    for panel = 1:2
        subplot(1,2,panel); hold on; grid on; axis equal;
        for e = 1:nElem
            xy = shape{e}.xy/1000;    % Display geometry in metres.
            h0 = plot(xy(1,:),xy(2,:),'k--','LineWidth',1.4);
            if panel == 1
                plot(xy(1,:),xy(2,:),'Color',[0.12 0.4 0.6],'LineWidth',3);
                midpoint = mean(xy,2);
                text(midpoint(1)+0.06,midpoint(2)+0.08,sprintf('e%d',e), ...
                    'Color',[0.15 0.35 0.6],'FontWeight','bold');
            else
                xyd = (shape{e}.xy+scale*shape{e}.displacement)/1000;
                h1 = plot(xyd(1,:),xyd(2,:),'-','Color',[0.85 0.2 0.12],'LineWidth',2);
            end
        end
        plot(nodes(:,1)/1000,nodes(:,2)/1000,'ko','MarkerFaceColor','w');
        for n = supportNodes.'
            if ismember(dofMap(n,3),fixedDofs)
                marker = 's'; label = 'Fixed';
            else
                marker = '^'; label = 'Pin';
            end
            plot(nodes(n,1)/1000,nodes(n,2)/1000,marker, ...
                'Color','k','MarkerSize',10,'LineWidth',1.5);
            text(nodes(n,1)/1000,nodes(n,2)/1000-0.13,label, ...
                'HorizontalAlignment','center','VerticalAlignment','top');
        end
        if panel == 1
            for n = 1:size(nodes,1)
                text(nodes(n,1)/1000+0.08,nodes(n,2)/1000+0.08, ...
                    sprintf('Node %d',n),'FontWeight','bold');
            end
            title([caption ' | model and loading']);
        else
            dd = reshape(d(dofMap.'),3,size(nodes,1)).';
            plot((nodes(:,1)+scale*dd(:,1))/1000, ...
                (nodes(:,2)+scale*dd(:,2))/1000,'o','Color',[0.85 0.2 0.12]);
            title(sprintf('Deformed shape | magnification %.1f',scale));
            legend([h0 h1],{'Undeformed','Deformed'},'Location','best');
        end
        xlabel('x (m)'); ylabel('y (m)');
        pad = 0.25*span/1000;
        xlim([min(nodes(:,1))/1000-pad,max(nodes(:,1))/1000+pad]);
        ylim([min(nodes(:,2))/1000-pad,max(nodes(:,2))/1000+pad]);
        set(gca,'FontSize',11);
    end
end
