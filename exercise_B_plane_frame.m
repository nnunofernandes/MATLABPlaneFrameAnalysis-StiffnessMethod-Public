%% EXERCISE B -- PLANE FRAME BY THE DIRECT STIFFNESS METHOD
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
% Nodes follow the figure: 1 left base, 2 left top, 3 right top, 4 right base.
% Elements: 1 = left column (up), 2 = beam (right), 3 = right column (down).
% The nonsequential dofMap below preserves the figure's DOF NUMBERS.
E = 200e3;                 % N/mm^2 = 200 GPa
I = 300e6;                 % mm^4
A = 9.68e3;                % mm^2
nodes = [0 0;             % [x y], mm
         0 3000;
         2000 3000;
         2000 0];
elements = [1 2; 2 3; 3 4];
dofMap = [11 12 8;         % Node 1: ux uy theta
          1  2 3;         % Node 2
          4  5 6;         % Node 3
          9 10 7];        % Node 4
fixedDofs = [9 10 11 12];  % Only translations at the two pin supports.
prescribedValues = zeros(4,1);
% IMPORTANT: DOFs 7 and 8 (base rotations) remain FREE. Do not set them to
% zero: that would create fixed supports. No internal element release is
% needed at a base connected to just one member; its free rotation enforces
% zero end moment. Upper nodes share one rotation each (rigid joints).

Fpoint = zeros(12,1);      % No concentrated nodal loads.

% 6 kN/m = 6 N/mm acts globally to the RIGHT along element 1.
% Element 1 points UP: c = 0, s = 1; its local +y points LEFT.
% Thus [qx; qy] = [c s; -s c]*[6;0] = [0;-6] N/mm.
w = 6;                    % Positive magnitude of rightward load, N/mm.
qLocal = [0 -w;
          0  0;
          0  0];
fprintf('EXERCISE B: reactions at pin supports 1 and 4\n');
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

%% 6. Independent analytical check: force method with one redundant
% Release only the horizontal restraint at node 4. The resulting primary
% frame has a left pin and a right vertical roller, and is determinate.
% Choose X = horizontal reaction at node 4, positive right.
%
% Unit-load integration gives horizontal displacement delta0 from the UDL:
% delta0 = w/(E*I)*(5*H^4/24 + H^3*B/4).
% The flexibility at this released DOF, including beam axial deformation:
% flexibility = (2*H^3/3 + H^2*B)/(E*I) + B/(E*A).
% Compatibility requires delta0 + flexibility*X = 0.
% Both columns contribute bending flexibility; the beam contributes bending
% and axial flexibility. Column axial terms vanish for this unit load.
B = 2000; H = 3000;
delta0 = w*(5*H^4/24 + H^3*B/4)/(E*I);
flexibility = (2*H^3/3 + H^2*B)/(E*I) + B/(E*A);
X = -delta0/flexibility;
Rx1 = -w*H-X;
Rx4 = X;
Ry4 = w*H^2/(2*B);
Ry1 = -Ry4;
rAnalytical = [Rx1; Ry1; Rx4; Ry4];
rComputed = R([dofMap(1,1:2) dofMap(4,1:2)]);
assert(all(abs(rComputed-rAnalytical) < 1e-8*max(abs(rAnalytical),1)), ...
    'Independent force-method check failed.');
% Pins cannot transmit moments. Local end i of member 1 is node 1;
% local end j of member 3 is node 4.
assert(abs(endForcesLocal(1,3)) < 1e-5, 'Left pin has a nonzero end moment.');
assert(abs(endForcesLocal(3,6)) < 1e-5, 'Right pin has a nonzero end moment.');
fprintf('\nREQUESTED ANSWER -- PIN SUPPORT REACTIONS\n');
fprintf('Node 1: Rx = %+.8f kN, Ry = %+.8f kN, M = 0 kN*m\n',Rx1/1e3,Ry1/1e3);
fprintf('Node 4: Rx = %+.8f kN, Ry = %+.8f kN, M = 0 kN*m\n',Rx4/1e3,Ry4/1e3);
fprintf('Negative Rx = left; negative Ry = down. A negative Ry is uplift demand.\n');
fprintf('Force-method and stiffness-method reactions agree.\n');

%% 7. Plot undeformed and deformed shapes
if makePlots
    plotFrame(nodes,elements,dofMap,d,elementData,qLocal,E,A,I, ...
        supportNodes,fixedDofs,deformationScale,'Exercise B');
    subplot(1,2,1); hold on;
    yy = linspace(0.1,2.9,10);
    quiver(-0.45*ones(size(yy)),yy,0.45*ones(size(yy)),zeros(size(yy)), ...
        0,'r','LineWidth',1.3,'MaxHeadSize',0.4);
    text(-0.48,1.5,'6 kN/m','Color','r','Rotation',90, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
end
% Useful workspace variables: K, F, Kff, d, R, elementData,
% displacementTable, reactionTable, endForceTable, rAnalytical.

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
