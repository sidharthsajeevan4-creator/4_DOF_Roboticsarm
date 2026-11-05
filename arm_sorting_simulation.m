1;
% =====================================================================
%  4-DOF ROBOTIC ARM  -  CONVEYOR PICK-AND-SORT SIMULATION  (GNU Octave)
%
%  A red cube, blue can and green ball ride a conveyor. The arm picks
%  each one and drops it into the bin of its own colour.
%
%  Arm data from the project deck (mm):
%     L1 = 120 (base)   L2 = 180 (upper arm)   L3 = 160 (forearm)
%     L4 = 80  (wrist)  gripper = 100
%
%  HOW TO RUN:   cd to this file's folder in Octave, then type:
%                arm_sorting_simulation
%
%  NOTE: use the default Octave GUI (qt graphics). The old "gnuplot"
%  toolkit cannot draw the shaded 3-D parts.
% =====================================================================

% ------------------------- KINEMATICS --------------------------------
% Forward kinematics: points base, shoulder, elbow, wrist joint,
% end of wrist link L4, gripper tip.
function P = fk(th, L, Lg)
  a2 = th(2);  a3 = th(2) + th(3);  a4 = th(2) + th(3) + th(4);
  rz = zeros(6, 2);                       % [radial distance, height]
  rz(2,:) = [0 L(1)];
  rz(3,:) = rz(2,:) + L(2) * [cos(a2) sin(a2)];
  rz(4,:) = rz(3,:) + L(3) * [cos(a3) sin(a3)];
  rz(5,:) = rz(4,:) + L(4) * [cos(a4) sin(a4)];
  rz(6,:) = rz(5,:) + Lg   * [cos(a4) sin(a4)];
  P = [rz(:,1) * cos(th(1)), rz(:,1) * sin(th(1)), rz(:,2)];
end

% Inverse kinematics for the gripper tip [x y z], gripper pointing down
function th = ik(p, L, Lg)
  th1 = atan2(p(2), p(1));
  r   = hypot(p(1), p(2));
  wz  = p(3) + Lg + L(4);
  h   = wz - L(1);
  D   = hypot(r, h);
  if D > L(2) + L(3) || D < abs(L(2) - L(3))
    error('Target [%g %g %g] is out of reach', p(1), p(2), p(3));
  end
  c3  = (D^2 - L(2)^2 - L(3)^2) / (2 * L(2) * L(3));
  th3 = -acos(c3);                                      % elbow up
  th2 = atan2(h, r) - atan2(L(3) * sin(th3), L(2) + L(3) * cos(th3));
  th4 = -pi/2 - th2 - th3;                              % gripper points down
  th  = [th1 th2 th3 th4];
end

function s = ease(u)
  s = (1 - cos(pi * u)) / 2;
end

% Cylindrical interpolation of the tip (straight lifts, circular swings)
function [TH, G] = move_cyl(p0, p1, n, g0, g1, L, Lg)
  c0 = [atan2(p0(2), p0(1)), hypot(p0(1), p0(2)), p0(3)];
  c1 = [atan2(p1(2), p1(1)), hypot(p1(1), p1(2)), p1(3)];
  TH = zeros(n, 4);  G = zeros(n, 1);
  for k = 1:n
    s = ease((k - 1) / (n - 1));
    c = c0 + s * (c1 - c0);
    TH(k,:) = ik([c(2)*cos(c(1)), c(2)*sin(c(1)), c(3)], L, Lg);
    G(k) = g0 + s * (g1 - g0);
  end
end

function [TH, G] = move_joint(th0, th1, n, g)
  TH = zeros(n, 4);
  for k = 1:n
    s = ease((k - 1) / (n - 1));
    TH(k,:) = th0 + s * (th1 - th0);
  end
  G = g * ones(n, 1);
end

function [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, belt, label)
  n = size(a, 1);
  if isscalar(belt), belt = belt * ones(n, 1); end
  TH = [TH; a];  G = [G; b];  BELT = [BELT; belt(:)];
  PH = [PH; repmat({label}, n, 1)];
end

% ------------------------- 3-D MESH TOOLS ----------------------------
% A "mesh" is a struct: V = vertices, F = quad faces, C = colour per face
function M = M_new(V, F, col)
  M.V = V;  M.F = F;  M.C = repmat(col, size(F, 1), 1);
end

function M = M_join(A, B)
  M.V = [A.V; B.V];
  M.F = [A.F; B.F + size(A.V, 1)];
  M.C = [A.C; B.C];
end

function M = M_box(x0, x1, y0, y1, z0, z1, col)
  V = [x0 y0 z0; x1 y0 z0; x1 y1 z0; x0 y1 z0; ...
       x0 y0 z1; x1 y0 z1; x1 y1 z1; x0 y1 z1];
  F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
  M = M_new(V, F, col);
end

% Surface of revolution around the z axis from a profile [radius height]
function M = M_revolve(prof, ns, col)
  m = size(prof, 1);
  ph = 2 * pi * (0:ns-1) / ns;
  V = zeros(m * ns, 3);
  for i = 1:m
    for j = 1:ns
      V((i-1)*ns + j, :) = [prof(i,1)*cos(ph(j)), prof(i,1)*sin(ph(j)), prof(i,2)];
    end
  end
  F = zeros((m-1) * ns, 4);  c = 0;
  for i = 1:m-1
    for j = 1:ns
      j2 = mod(j, ns) + 1;
      c = c + 1;
      F(c,:) = [(i-1)*ns + j, (i-1)*ns + j2, i*ns + j2, i*ns + j];
    end
  end
  M = M_new(V, F, col);
end

function M = M_cyl(r, z0, z1, col, ns)
  M = M_revolve([0 z0; r z0; r z1; 0 z1], ns, col);
end

function M = M_sphere(R, col)
  a = linspace(0, pi, 11)';
  M = M_revolve([R * sin(a), -R * cos(a)], 20, col);
end

% Simple directional shading: brighter faces point towards the light
function C = shade(V, F, col)
  A  = V(F(:,3), :) - V(F(:,1), :);
  Bv = V(F(:,4), :) - V(F(:,2), :);
  N  = [A(:,2).*Bv(:,3) - A(:,3).*Bv(:,2), ...
        A(:,3).*Bv(:,1) - A(:,1).*Bv(:,3), ...
        A(:,1).*Bv(:,2) - A(:,2).*Bv(:,1)];
  nn = sqrt(sum(N.^2, 2));
  nn(nn < 1e-9) = 1;
  N  = N ./ nn;
  Ld = [-0.35 -0.55 0.75];  Ld = Ld / norm(Ld);
  k  = 0.38 + 0.62 * abs(N * Ld');
  C  = min(1, k .* col);
end

% Rotation matrix whose 3rd column is direction d
function R = frame_from_dir(d)
  d = d(:)' / norm(d);
  if abs(d(3)) > 0.999
    xp = [1 0 0];
  else
    xp = cross([0 0 1], d);  xp = xp / norm(xp);
  end
  yp = cross(d, xp);
  R = [xp(:) yp(:) d(:)];
end

function V = xf(V0, R, t)
  V = V0 * R' + t;
end

function h = draw_mesh(M, V)
  h = patch('Vertices', V, 'Faces', M.F, 'FaceVertexCData', shade(V, M.F, M.C), ...
            'FaceColor', 'flat', 'EdgeColor', 'none');
end

function set_mesh(h, M, V)
  set(h, 'Vertices', V, 'FaceVertexCData', shade(V, M.F, M.C));
end

function [X, Y, Z] = belt_stripes(belt)
  xs = -575 + mod((0:13) * 50 + belt, 685);
  X = [xs; xs; nan(1, numel(xs))];  X = X(:)';
  Y = repmat([-277 -223 NaN], 1, numel(xs));
  Z = 50.6 * ones(size(Y));
end

% =====================================================================
%  MAIN PROGRAM
% =====================================================================
clc; close all;
if exist('OCTAVE_VERSION', 'builtin')
  try, graphics_toolkit('qt'); catch, end
end

% ---- Arm and task settings (mm) ----
L  = [120 180 160 80];  Lg = 100;
home    = [0 90 -90 -90] * pi / 180;
gOpen   = 72;
beltTop = 50;   binFloor = 8;   hoverUp = 60;   binHoverZ = 150;
pickXY  = [0 -250];

% item types: 1 = red cube, 2 = blue can, 3 = green ball
types   = [1 2 3 1];
nItems  = numel(types);
Hitem   = [36 44 44];          % item heights
Gclosed = [36 40 44];          % gripper width when holding each type
binAng  = [-45 0 45];          % bin angle around the arm base (deg)
binC    = 250 * [cosd(binAng') sind(binAng')];
colItem = [0.88 0.18 0.14; 0.15 0.42 0.90; 0.20 0.72 0.30];
binName = {'RED CUBES', 'BLUE CANS', 'GREEN BALLS'};

B  = 200 + 100 * (0:nItems-1);     % conveyor position when item i is at the pick point
x0 = -B;                           % where each item starts on the belt

% ---- Build the complete motion ----
TH = zeros(0, 4);  G = zeros(0, 1);  BELT = zeros(0, 1);  PH = cell(0, 1);
holdStart = zeros(1, nItems);  holdEnd = zeros(1, nItems);
relXY = zeros(nItems, 2);  restZ = zeros(1, nItems);  zRel = zeros(1, nItems);

s = ease(linspace(0, 1, 50)');
[TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, repmat(home, 50, 1), gOpen * ones(50, 1), ...
                               B(1) * s, 'Conveyor feeding items');

for ii = 1:nItems
  t = types(ii);  h = Hitem(t);
  c = sum(types(1:ii-1) == t);                       % items already in this bin
  tv = [-sind(binAng(t)) cosd(binAng(t))];
  placeXY  = binC(t,:) + (-1)^c * 22 * tv;
  restZ(ii) = binFloor + h / 2;
  zRel(ii)  = restZ(ii) + 25;
  relXY(ii,:) = placeXY;

  pPick    = [pickXY, beltTop + h/2];
  pPickUp  = [pickXY, beltTop + h/2 + hoverUp];
  pPlaceUp = [placeXY, binHoverZ];
  pPlace   = [placeXY, zRel(ii)];
  gC = Gclosed(t);
  nm = binName{t};

  [a, b] = move_joint(home, ik(pPickUp, L, Lg), 45, gOpen);
  [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, B(ii), 'Moving to pick position');
  [a, b] = move_cyl(pPickUp, pPick, 22, gOpen, gOpen, L, Lg);
  [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, B(ii), 'Descending to item');
  [a, b] = move_cyl(pPick, pPick, 12, gOpen, gC, L, Lg);
  [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, B(ii), 'Closing gripper');
  holdStart(ii) = size(TH, 1) + 1;
  [a, b] = move_cyl(pPick, pPickUp, 22, gC, gC, L, Lg);
  [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, B(ii), 'Lifting item');
  [a, b] = move_cyl(pPickUp, pPlaceUp, 45, gC, gC, L, Lg);
  [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, B(ii), ['Swinging to ' nm ' bin']);
  [a, b] = move_cyl(pPlaceUp, pPlace, 22, gC, gC, L, Lg);
  [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, B(ii), 'Lowering into bin');
  holdEnd(ii) = size(TH, 1);
  [a, b] = move_cyl(pPlace, pPlace, 12, gC, gOpen, L, Lg);
  [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, B(ii), 'Releasing item');
  [a, b] = move_cyl(pPlace, pPlaceUp, 18, gOpen, gOpen, L, Lg);
  [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, B(ii), 'Lifting away');

  if ii < nItems                                    % belt advances while arm returns
    bl = B(ii) + (B(ii+1) - B(ii)) * ease(linspace(0, 1, 40)');
    lab = 'Returning home / next item arriving';
  else
    bl = B(ii);  lab = 'Returning home';
  end
  [a, b] = move_joint(ik(pPlaceUp, L, Lg), home, 40, gOpen);
  [TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, a, b, bl, lab);
end
[TH, G, BELT, PH] = add_frames(TH, G, BELT, PH, repmat(home, 40, 1), gOpen * ones(40, 1), ...
                               B(nItems), 'Sorting complete');
N = size(TH, 1);

TIP = zeros(N, 3);
for k = 1:N
  P = fk(TH(k,:), L, Lg);  TIP(k,:) = P(6,:);
end
printf('Planned %d frames for %d items. Opening 3-D window...\n', N, nItems);

% ---- Colours ----
cWhite  = [0.93 0.94 0.96];  cOrange = [0.96 0.46 0.10];
cDark   = [0.17 0.19 0.23];  cChrome = [0.72 0.74 0.78];
binCol  = [0.78 0.22 0.18; 0.18 0.38 0.80; 0.18 0.60 0.28];

% ---- Figure and axes ----
fig = figure('Name', '4-DOF Arm - Pick and Sort', 'NumberTitle', 'off', ...
             'Color', [0.05 0.06 0.09], 'Position', [60 60 1150 720]);
ax = axes('Position', [0 0 1 1]);
hold on;
axis([-600 400 -400 400 0 520]);
daspect([1 1 1]);
axis off;
view(-30, 22);
try, set(ax, 'Projection', 'perspective'); catch, end

% ---- Floor (checker tiles) ----
fV = [];  fF = [];  fC = [];  c = 0;
for xi = -600:100:300
  for yi = -400:100:300
    v0 = size(fV, 1);
    fV = [fV; xi yi 0; xi+100 yi 0; xi+100 yi+100 0; xi yi+100 0];
    fF = [fF; v0+1 v0+2 v0+3 v0+4];
    if mod((xi + 600) / 100 + (yi + 400) / 100, 2) == 0
      fC = [fC; 0.21 0.23 0.28];
    else
      fC = [fC; 0.27 0.29 0.35];
    end
  end
end
patch('Vertices', fV, 'Faces', fF, 'FaceVertexCData', fC, 'FaceColor', 'flat', 'EdgeColor', 'none');

% ---- Walls with a vertical gradient ----
wb = [0.30 0.34 0.42];  wt = [0.09 0.11 0.16];
patch('Vertices', [-600 400 0; 400 400 0; 400 400 520; -600 400 520], 'Faces', [1 2 3 4], ...
      'FaceVertexCData', [wb; wb; wt; wt], 'FaceColor', 'interp', 'EdgeColor', 'none');
patch('Vertices', [400 -400 0; 400 400 0; 400 400 520; 400 -400 520], 'Faces', [1 2 3 4], ...
      'FaceVertexCData', [wb; wb; wt; wt], 'FaceColor', 'interp', 'EdgeColor', 'none');
patch('Vertices', [-600 399 130; 400 399 130; 400 399 142; -600 399 142], 'Faces', [1 2 3 4], ...
      'FaceVertexCData', repmat([0.10 0.65 0.80], 4, 1), 'FaceColor', 'flat', 'EdgeColor', 'none');
patch('Vertices', [399 -400 130; 399 400 130; 399 400 142; 399 -400 142], 'Faces', [1 2 3 4], ...
      'FaceVertexCData', repmat([0.10 0.65 0.80], 4, 1), 'FaceColor', 'flat', 'EdgeColor', 'none');
for px = -450:150:-0                                  % wall panels
  plot3([px px+110 px+110 px px], [398 398 398 398 398], [190 190 330 330 190], ...
        '-', 'Color', [0.45 0.50 0.60], 'LineWidth', 1);
end
text(-420, 396, 420, 'ROBOTICS LAB  |  4-DOF PICK & SORT CELL', 'Color', [0.85 0.92 1.0], ...
     'FontSize', 15, 'FontWeight', 'bold');

% ---- Static equipment: pedestal, conveyor, bins ----
S = M_cyl(105, 0.3, 0.6, [0.10 0.10 0.12], 40);                    % floor shadow
S = M_join(S, M_revolve([175 0.7; 183 0.7], 60, [0.95 0.78 0.10])); % safety ring
S = M_join(S, M_cyl(95, 0, 14, cDark, 32));                        % base plate
S = M_join(S, M_cyl(62, 14, 70, [0.26 0.28 0.33], 32));            % column

S = M_join(S, M_box(-580, 110, -280, -220, 38, 50, [0.10 0.10 0.12]));   % belt
S = M_join(S, M_box(-580, 110, -278, -222, 28, 38, cDark));              % under-frame
S = M_join(S, M_box(-585, 115, -292, -280, 30, 58, cChrome));            % rails
S = M_join(S, M_box(-585, 115, -220, -208, 30, 58, cChrome));
S = M_join(S, M_box(115, 160, -275, -225, 20, 62, cOrange));             % drive motor
for lx = [-570 -240 95]
  for ly = [-292 -208]
    S = M_join(S, M_box(lx, lx+12, ly, ly+12, 0, 30, [0.45 0.47 0.52]));
  end
end
Rr = frame_from_dir([0 1 0]);
for rx = [-580 110]                                                      % end rollers
  Mr = M_cyl(14, -38, 38, [0.55 0.57 0.62], 16);
  Mr.V = xf(Mr.V, Rr, [rx -250 44]);
  S = M_join(S, Mr);
end
for bi = 1:3                                                             % bins
  cx = binC(bi,1);  cy = binC(bi,2);  bc = binCol(bi,:);
  S = M_join(S, M_box(cx-60, cx+60, cy-60, cy+60, 0, binFloor, bc * 0.7));
  S = M_join(S, M_box(cx-60, cx+60, cy-60, cy-54, 0, 60, bc));
  S = M_join(S, M_box(cx-60, cx+60, cy+54, cy+60, 0, 60, bc));
  S = M_join(S, M_box(cx-60, cx-54, cy-54, cy+54, 0, 60, bc));
  S = M_join(S, M_box(cx+54, cx+60, cy-54, cy+54, 0, 60, bc));
  text(cx, cy, 78, binName{bi}, 'Color', [1 1 1], 'FontSize', 9, ...
       'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end
draw_mesh(S, S.V);
text(-300, -250, 90, 'CONVEYOR', 'Color', [0.9 0.9 0.95], 'FontSize', 9, 'HorizontalAlignment', 'center');

% ---- Moving belt stripes ----
[sx, sy, sz] = belt_stripes(BELT(1));
hStripe = plot3(sx, sy, sz, '-', 'Color', [0.38 0.40 0.45], 'LineWidth', 2);

% ---- Robot arm parts (local meshes) ----
Mturret = M_join(M_cyl(56, 70, 110, cOrange, 32), ...
                 M_join(M_box(-34, 34, 30, 42, 80, 140, cOrange), ...
                 M_join(M_box(-34, 34, -42, -30, 80, 140, cOrange), ...
                        M_box(-64, -42, -22, 22, 78, 112, cDark))));
Mup    = M_cyl(21, 0, 1, cWhite, 20);
Mfore  = M_cyl(16, 0, 1, cWhite, 20);
Mwlink = M_cyl(12, 0, 1, cDark, 16);
Msh    = M_cyl(30, -30, 30, cDark, 24);
Mel    = M_cyl(26, -30, 30, cOrange, 24);
Mwr    = M_cyl(20, -25, 25, cDark, 20);
Mpalm  = M_join(M_box(-16, 16, -50, 50, 0, 16, cDark), M_cyl(12, -6, 0, cChrome, 12));
MfP    = M_box(-10, 10, 0, 9, 10, Lg + 12, cChrome);
MfM    = M_box(-10, 10, -9, 0, 10, Lg + 12, cChrome);

P  = fk(TH(1,:), L, Lg);
hTurret = draw_mesh(Mturret, Mturret.V);
hUp   = draw_mesh(Mup, Mup.V);    hFore = draw_mesh(Mfore, Mfore.V);
hWl   = draw_mesh(Mwlink, Mwlink.V);
hSh   = draw_mesh(Msh, Msh.V);    hEl = draw_mesh(Mel, Mel.V);   hWr = draw_mesh(Mwr, Mwr.V);
hPalm = draw_mesh(Mpalm, Mpalm.V);
hFP   = draw_mesh(MfP, MfP.V);    hFM = draw_mesh(MfM, MfM.V);

% ---- Items ----
Mitem = cell(1, 3);
Mitem{1} = M_box(-18, 18, -18, 18, -18, 18, colItem(1,:));
Mitem{2} = M_cyl(20, -22, 22, colItem(2,:), 24);
Mitem{3} = M_sphere(22, colItem(3,:));
hItem = zeros(1, nItems);
for ii = 1:nItems
  hItem(ii) = draw_mesh(Mitem{types(ii)}, Mitem{types(ii)}.V);
end

hTrail = plot3(TIP(1,1), TIP(1,2), TIP(1,3), '-', 'Color', [0.2 0.9 1.0], 'LineWidth', 1.2);

% ---- Heads-up display ----
hHud = axes('Position', [0 0.93 1 0.07], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1], 'Color', 'none');
hTxt = text(0.02, 0.5, '', 'Parent', hHud, 'Color', [1 1 1], 'FontSize', 12, 'FontWeight', 'bold');
set(fig, 'CurrentAxes', ax);

% ---- Animate ----
frameDelay = 0.012;
for k = 1:N
  th = TH(k,:);  P = fk(th, L, Lg);  g = G(k);
  tvec = [-sin(th(1)) cos(th(1)) 0];                  % sideways direction of the arm plane

  % turret rotates with the base joint
  Rz = [cos(th(1)) -sin(th(1)) 0; sin(th(1)) cos(th(1)) 0; 0 0 1];
  set_mesh(hTurret, Mturret, xf(Mturret.V, Rz, [0 0 0]));

  % links (cylinders stretched between joints)
  d = P(3,:) - P(2,:);  Vs = Mup.V;     Vs(:,3) = Vs(:,3) * norm(d);
  set_mesh(hUp, Mup, xf(Vs, frame_from_dir(d), P(2,:)));
  d = P(4,:) - P(3,:);  Vs = Mfore.V;   Vs(:,3) = Vs(:,3) * norm(d);
  set_mesh(hFore, Mfore, xf(Vs, frame_from_dir(d), P(3,:)));
  d = P(5,:) - P(4,:);  Vs = Mwlink.V;  Vs(:,3) = Vs(:,3) * norm(d);
  set_mesh(hWl, Mwlink, xf(Vs, frame_from_dir(d), P(4,:)));

  % joint housings (axis = sideways direction)
  Rt = frame_from_dir(tvec);
  set_mesh(hSh, Msh, xf(Msh.V, Rt, P(2,:)));
  set_mesh(hEl, Mel, xf(Mel.V, Rt, P(3,:)));
  set_mesh(hWr, Mwr, xf(Mwr.V, Rt, P(4,:)));

  % gripper: z' = pointing direction, y' = sideways
  zp = (P(6,:) - P(5,:)) / Lg;
  yp = tvec;
  xp = cross(yp, zp);
  Rg = [xp(:) yp(:) zp(:)];
  set_mesh(hPalm, Mpalm, xf(Mpalm.V, Rg, P(5,:)));
  V = MfP.V;  V(:,2) = V(:,2) + g / 2;  set_mesh(hFP, MfP, xf(V, Rg, P(5,:)));
  V = MfM.V;  V(:,2) = V(:,2) - g / 2;  set_mesh(hFM, MfM, xf(V, Rg, P(5,:)));

  % items: on the belt, in the gripper, or sorted in a bin
  for ii = 1:nItems
    h = Hitem(types(ii));
    if k < holdStart(ii)
      pos = [x0(ii) + BELT(k), pickXY(2), beltTop + h/2];
    elseif k <= holdEnd(ii)
      pos = P(6,:);
    else
      pos = [relXY(ii,:), max(restZ(ii), zRel(ii) - 8 * (k - holdEnd(ii)))];
    end
    Mi = Mitem{types(ii)};
    set_mesh(hItem(ii), Mi, Mi.V + pos);
  end

  % conveyor stripes and tool path
  [sx, sy, sz] = belt_stripes(BELT(k));
  set(hStripe, 'XData', sx, 'YData', sy, 'ZData', sz);
  idx = max(1, k - 90):k;
  set(hTrail, 'XData', TIP(idx,1), 'YData', TIP(idx,2), 'ZData', TIP(idx,3));

  done = holdEnd < k;
  set(hTxt, 'String', sprintf('%s    |    Sorted %d/%d    |    th1 %4.0f   th2 %4.0f   th3 %4.0f   th4 %4.0f deg', ...
      PH{k}, sum(done), nItems, th(1)*180/pi, th(2)*180/pi, th(3)*180/pi, th(4)*180/pi));

  view(-30 + 6 * sin(2 * pi * k / N), 22);
  drawnow;
  pause(frameDelay);
end

% ---- Joint angle history ----
figure('Name', 'Joint angles', 'NumberTitle', 'off');
plot(1:N, TH * 180 / pi, 'LineWidth', 1.5);
grid on;  xlabel('Frame');  ylabel('Joint angle (deg)');
legend('th1 base', 'th2 shoulder', 'th3 elbow', 'th4 wrist', 'Location', 'best');
title('Joint angles during the sorting cycle');
disp('Done.');
