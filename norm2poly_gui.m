function norm2poly_gui(Z, fname)
%NORM2POLY_GUI  Interactive, tabbed GUI version of norm2poly.m
%
%   norm2poly_gui                  open an empty GUI
%   norm2poly_gui('Mn')            open the GUI with Mn preselected
%   norm2poly_gui('Mn','scan.txt') open the GUI and load scan.txt
%   norm2poly_gui(25,'scan.txt')   same, using the atomic number
%
%   Each tab is a completely independent dataset: its own file, element,
%   edge, fit windows, polynomial order and result.  Nothing is shared
%   between tabs, so a fit in one cannot disturb another.
%
%   The three buttons at the top of every tab (also on the "Datasets"
%   menu) manage them:
%     + New       an empty tab, everything at its defaults
%     Duplicate   a new tab carrying this tab's element, edge, edge
%                 energy, fit windows, polynomial order, unit and column
%                 choice, but no data - load the next scan into it and
%                 the settings are already right
%     Close       close this tab; closing the last one leaves a fresh
%                 empty tab behind
%
%   A duplicated tab keeps the copied fit windows for the first file you
%   load into it, provided they fall inside that scan; after that it
%   behaves like any other tab.
%
%   The fit is the one from norm2poly.m.  Over the pre-edge and post-edge
%   windows the raw signal is modelled as
%
%       y = c0*mu_McMaster + c1*P_0 + c2*P_1 + ... + cn*P_n
%
%   with P_j the Legendre polynomials on [-1,1] and mu_McMaster the
%   McMaster edge step after its extrapolated pre-edge baseline has been
%   removed.  The fitted polynomial is subtracted over the whole scan and
%   the result is put on the McMaster scale, [cm^2/g].
%
%   This file is self-contained.  McMaster.m, zymbol.m, unify.m,
%   legpoly.m and the element.mat tables are all embedded as local
%   functions at the bottom, so nothing else needs to be on the path.
%   The local copies are used in preference to any on the path.
%
%   Working in a tab
%     1. load a file - the raw scan is plotted straight away
%     2. give the element and pick the edge; the edge energy comes from
%        the McMaster tables but you can overwrite it
%     3. set the four window energies: type them, drag the dashed lines
%        on the top plot, or use the four sliders underneath it.  The
%        pre-edge window is green, the post-edge window blue.
%     4. the fit follows every change while "auto update" is ticked;
%        untick it and use "Run normalization" instead if you would
%        rather fit on demand
%
%   Each graph has its own zoom: the + - and A buttons above it, the
%   mouse wheel over it, or drag out a box on it to zoom into exactly
%   that region.  A restores the full view.  A zoom is kept when the fit
%   is recomputed, and is reset when a new file is loaded.
%
%   The two graphs share the vertical space; drag the grey bar between
%   them to give one more room.  "Graphs -> Pop out" copies either one
%   into an ordinary figure window for zooming, printing or saving.
%
%   "To workspace" assigns a struct to the base workspace holding the
%   energy, the normalized signal, the McMaster reference, the fitted
%   background, the coefficients and the fit windows - everything the
%   plots are drawn from, for your own analysis at the command line.
%
%   Based on norm2poly.m, Tsu-Chien Weng, Penner-Hahn research group.

if nargin < 1, Z = []; end
if nargin < 2, fname = ''; end

% Built hidden: every uicontrol added to a visible figure forces a
% redraw, and there are about fifty of them, which is what made start-up
% slow.  Shown in one go at the end instead.
%
% ToolBar 'none' drops the icon strip, whose buttons each cost time to
% build and all duplicate an entry in File / Edit / View / Insert /
% Tools.  Turn it back on for a session with View -> Figure Toolbar if
% you want it.  DockControls 'off' saves the docking arrow as well.
fig = figure('Name','norm2poly  -  XAS background removal', ...
    'NumberTitle','off','Units','normalized','Visible','off', ...
    'ToolBar','none','DockControls','off', ...
    'Position',[0.08 0.10 0.80 0.80],'Color',[0.94 0.94 0.94]);

setappdata(fig,'element',n2p_elementTable());   % shared, read only

m = uimenu(fig,'Label','Datasets');
uimenu(m,'Label','New tab','Accelerator','T','Callback',@(s,e)menuNewTab(fig));
uimenu(m,'Label','Duplicate tab','Accelerator','D','Callback',@(s,e)menuDuplicateTab(fig));
uimenu(m,'Label','Rename tab...','Callback',@(s,e)menuRenameTab(fig));
uimenu(m,'Label','Close tab','Separator','on','Callback',@(s,e)menuCloseTab(fig));

g = uimenu(fig,'Label','Graphs');
uimenu(g,'Label','Pop out top graph','Callback',@(s,e)popOut(fig,1));
uimenu(g,'Label','Pop out bottom graph','Callback',@(s,e)popOut(fig,2));
uimenu(g,'Label','Reset graph split','Separator','on','Callback',@(s,e)resetSplit(fig));

tg = uitabgroup('Parent',fig,'Units','normalized','Position',[0 0 1 1]);
setappdata(fig,'tabgroup',tg);
tb = newTab(fig);
armHover(fig);
set(fig,'WindowScrollWheelFcn',@(a,b)onScroll(fig,b));

if ~isempty(Z)
    s = st(tb);
    if ischar(Z), set(s.edZ,'String',Z); else set(s.edZ,'String',num2str(Z)); end
    onElement(tb);
end
if ~isempty(fname), loadFile(tb,fname); end

set(fig,'Visible','on');
drawnow;
end

% ======================================================================
%                          TABS AND STATE
% ======================================================================
function s = st(tb)
s = getappdata(tb,'state');
end

function stSet(tb,s)
setappdata(tb,'state',s);
end

function el = elements(s)
el = getappdata(s.fig,'element');
end

% ----------------------------------------------------------------------
function tb = newTab(fig,name)
tg = getappdata(fig,'tabgroup');
if nargin < 2
    name = sprintf('dataset %d',numel(get(tg,'Children'))+1);
end

tb = uitab('Parent',tg,'Title',name,'BackgroundColor',[0.94 0.94 0.94]);
s  = buildTab(fig,tb);

s.data = [];  s.file = '';  s.path = pwd;
s.x = [];  s.y = [];  s.Z = [];  s.result = [];
s.stale = true;  s.drag = 0;
s.wantCols = [1 2];      % preferred energy / signal columns
s.zoomed = [false false];% true once a graph has been zoomed by hand
s.boxK = 0;              % graph currently being box-zoomed, 0 = none
s.boxAnchor = [NaN NaN];
s.boxH = gobjects(1,0);
s.curveH = {gobjects(1,0), gobjects(1,0)};   % curves the cursor can track
s.curveName = {{},{}};
s.keepRange = false;     % set on a duplicate: keep the copied windows
s.line = gobjects(1,4);  s.patch = gobjects(1,2);
setappdata(tb,'lock',0);
stSet(tb,s);

set(tg,'SelectedTab',tb);
end

% ----------------------------------------------------------------------
function menuNewTab(fig)
newTab(fig);
end

function menuDuplicateTab(fig)
tg = getappdata(fig,'tabgroup');
tb = get(tg,'SelectedTab');
if ~isempty(tb), duplicateTab(tb); end
end

% ----------------------------------------------------------------------
function tb2 = duplicateTab(tb)
%DUPLICATETAB  A new tab set up exactly like this one, but with no data.
%   Everything that describes how to process a scan is carried over; the
%   scan itself is not, so the copy is ready for the next file.
s = st(tb);

title0 = get(tb,'Title');
tb2 = newTab(s.fig,[title0 ' (2)']);
s2  = st(tb2);

set(s2.popUnit, 'Value', get(s.popUnit,'Value'));
set(s2.edZ,     'String',get(s.edZ,'String'));
set(s2.popEdge, 'Value', get(s.popEdge,'Value'));
set(s2.edE0,    'String',get(s.edE0,'String'));
set(s2.popOrder,'Value', get(s.popOrder,'Value'));
set(s2.chkAuto, 'Value', get(s.chkAuto,'Value'));

r = getRange(s);
if all(isfinite(r))
    set(s2.edPre1,'String',get(s.edPre1,'String'));
    set(s2.edPre2,'String',get(s.edPre2,'String'));
    set(s2.edPost1,'String',get(s.edPost1,'String'));
    set(s2.edPost2,'String',get(s.edPost2,'String'));
    s2.keepRange = true;
end

s2.Z        = s.Z;
s2.path     = s.path;
s2.wantCols = [get(s.popEcol,'Value') get(s.popYcol,'Value')];
stSet(tb2,s2);

setStatus(s2,{'Settings copied from',['"' title0 '".'],'Load a data file.'});
end

function menuCloseTab(fig)
tg = getappdata(fig,'tabgroup');
tb = get(tg,'SelectedTab');
if ~isempty(tb), closeTab(tb); end
end

function closeTab(tb)
if ~ishghandle(tb), return; end
s   = st(tb);
fig = s.fig;
if ishghandle(fig)      % a drag in this tab must not outlive it
    armHover(fig);
end
tg = getappdata(fig,'tabgroup');
delete(tb);
if isempty(get(tg,'Children')), newTab(fig); end
end

function menuRenameTab(fig)
tg = getappdata(fig,'tabgroup');
tb = get(tg,'SelectedTab');
if isempty(tb), return; end
answer = inputdlg('Tab name:','Rename tab',1,{get(tb,'Title')});
if isempty(answer), return; end
if ~isempty(strtrim(answer{1})), set(tb,'Title',strtrim(answer{1})); end
end

% ======================================================================
%                             INTERFACE
% ======================================================================
function s = buildTab(fig,tb)

bg   = [0.94 0.94 0.94];
GRN  = [0.10 0.55 0.25];       % pre-edge colour
BLU  = [0.10 0.35 0.85];       % post-edge colour
hRow = 0.040;  hLbl = 0.028;  gap = 0.006;  sec = 0.012;
x0   = 0.05;  w = 0.90;

s.fig = fig;  s.tab = tb;
s.colPre = GRN;  s.colPost = BLU;

pnl = uipanel('Parent',tb,'Units','normalized','Visible','off', ...
    'Position',[0.005 0.005 0.245 0.990],'BackgroundColor',bg);
s.pnl = pnl;                 % made visible again at the end of buildTab

top = 0.985;

% ---------------- tab management ----------------
uicontrol(pnl,'Style','pushbutton','Units','normalized', ...
    'Position',[x0 top-0.034 0.27 0.034],'String','+ New', ...
    'TooltipString','open an empty tab', ...
    'Callback',@(a,b)newTab(fig));
uicontrol(pnl,'Style','pushbutton','Units','normalized', ...
    'Position',[x0+0.30 top-0.034 0.34 0.034],'String','Duplicate', ...
    'TooltipString','new tab with these settings, ready for the next file', ...
    'Callback',@(a,b)duplicateTab(tb));
uicontrol(pnl,'Style','pushbutton','Units','normalized', ...
    'Position',[x0+0.67 top-0.034 0.23 0.034],'String','Close', ...
    'TooltipString','close this tab', ...
    'Callback',@(a,b)closeTab(tb));
top = top-0.034-gap;

% ---------------- 1. data file ----------------
head(pnl,[x0 top-hLbl w hLbl],'1.  Data file',bg);        top = top-hLbl-gap;

uicontrol(pnl,'Style','pushbutton','Units','normalized', ...
    'Position',[x0 top-hRow w hRow],'String','Load data file...', ...
    'Callback',@(a,b)loadFile(tb));                       top = top-hRow-gap;

s.txFile = uicontrol(pnl,'Style','text','Units','normalized', ...
    'Position',[x0 top-hLbl w hLbl],'String','(no file loaded)', ...
    'HorizontalAlignment','left','ForegroundColor',[0.35 0.35 0.35], ...
    'BackgroundColor',bg);                                top = top-hLbl-gap;

lbl(pnl,[x0 top-hRow+0.004 0.28 hLbl],'energy col',bg);
s.popEcol = uicontrol(pnl,'Style','popupmenu','Units','normalized', ...
    'Position',[x0+0.28 top-hRow 0.14 hRow],'String',{'1'}, ...
    'BackgroundColor','w','Callback',@(a,b)onColumns(tb));
lbl(pnl,[x0+0.47 top-hRow+0.004 0.26 hLbl],'signal col',bg);
s.popYcol = uicontrol(pnl,'Style','popupmenu','Units','normalized', ...
    'Position',[x0+0.76 top-hRow 0.14 hRow],'String',{'2'}, ...
    'BackgroundColor','w','Callback',@(a,b)onColumns(tb));
top = top-hRow-gap;

lbl(pnl,[x0 top-hRow+0.004 0.34 hLbl],'energy unit',bg);
s.popUnit = uicontrol(pnl,'Style','popupmenu','Units','normalized', ...
    'Position',[x0+0.36 top-hRow 0.54 hRow],'String',{'eV','keV','auto'}, ...
    'Value',1,'BackgroundColor','w','Callback',@(a,b)onColumns(tb));
top = top-hRow-gap-sec;

% ---------------- 2. element and edge ----------------
head(pnl,[x0 top-hLbl w hLbl],'2.  Element and edge',bg); top = top-hLbl-gap;

lbl(pnl,[x0 top-hRow+0.004 0.20 hLbl],'element',bg);
s.edZ = uicontrol(pnl,'Style','edit','Units','normalized', ...
    'Position',[x0+0.21 top-hRow 0.20 hRow],'String','', ...
    'BackgroundColor','w','Callback',@(a,b)onElement(tb));
lbl(pnl,[x0+0.46 top-hRow+0.004 0.14 hLbl],'edge',bg);
s.popEdge = uicontrol(pnl,'Style','popupmenu','Units','normalized', ...
    'Position',[x0+0.60 top-hRow 0.30 hRow],'String',{'K','L1','L2','L3'}, ...
    'BackgroundColor','w','Callback',@(a,b)onEdge(tb));
top = top-hRow-gap;

lbl(pnl,[x0 top-hRow+0.004 0.52 hLbl],'edge energy [eV]',bg);
s.edE0 = uicontrol(pnl,'Style','edit','Units','normalized', ...
    'Position',[x0+0.54 top-hRow 0.36 hRow],'String','', ...
    'BackgroundColor','w','Callback',@(a,b)onE0(tb));
top = top-hRow-gap-sec;

% ---------------- 3. fit windows ----------------
head(pnl,[x0 top-hLbl w hLbl],'3.  Fit windows  [eV]',bg); top = top-hLbl-gap;

s.edPre1  = pairEdit(pnl,[x0 top-hRow w hRow],'pre-edge start', bg,GRN,@(a,b)onRange(tb));
top = top-hRow-gap;
s.edPre2  = pairEdit(pnl,[x0 top-hRow w hRow],'pre-edge end',   bg,GRN,@(a,b)onRange(tb));
top = top-hRow-gap;
s.edPost1 = pairEdit(pnl,[x0 top-hRow w hRow],'post-edge start',bg,BLU,@(a,b)onRange(tb));
top = top-hRow-gap;
s.edPost2 = pairEdit(pnl,[x0 top-hRow w hRow],'post-edge end',  bg,BLU,@(a,b)onRange(tb));
top = top-hRow-gap;

uicontrol(pnl,'Style','pushbutton','Units','normalized', ...
    'Position',[x0 top-hRow w hRow],'String','reset to  E0-20 / E0+80 eV', ...
    'Callback',@(a,b)defaultRange(tb));                   top = top-hRow-gap;

lbl(pnl,[x0 top-hRow+0.004 0.55 hLbl],'polynomial order',bg);
s.popOrder = uicontrol(pnl,'Style','popupmenu','Units','normalized', ...
    'Position',[x0+0.57 top-hRow 0.33 hRow], ...
    'String',{'0','1','2','3','4','5','6','7','8'},'Value',3, ...
    'BackgroundColor','w','Callback',@(a,b)afterParamChange(tb));
top = top-hRow-gap-sec;

% ---------------- 4. run ----------------
head(pnl,[x0 top-hLbl w hLbl],'4.  Normalize',bg);        top = top-hLbl-gap;

s.chkAuto = uicontrol(pnl,'Style','checkbox','Units','normalized', ...
    'Position',[x0 top-hRow w hRow],'String','auto update','Value',1, ...
    'BackgroundColor',bg,'Callback',@(a,b)onAuto(tb));    top = top-hRow-gap;

uicontrol(pnl,'Style','pushbutton','Units','normalized', ...
    'Position',[x0 top-hRow w hRow],'String','Run normalization', ...
    'FontWeight','bold','Callback',@(a,b)runFit(tb,true)); top = top-hRow-gap;

uicontrol(pnl,'Style','pushbutton','Units','normalized', ...
    'Position',[x0 top-hRow 0.43 hRow],'String','Save .nrm', ...
    'TooltipString','write energy, normalized and McMaster as ASCII', ...
    'Callback',@(a,b)saveNrm(tb));
uicontrol(pnl,'Style','pushbutton','Units','normalized', ...
    'Position',[x0+0.47 top-hRow 0.43 hRow],'String','To workspace', ...
    'TooltipString','put the whole fit in a struct in the base workspace', ...
    'Callback',@(a,b)toWorkspace(tb));                    top = top-hRow-gap;

% ---------------- status ----------------
s.txStat = uicontrol(pnl,'Style','text','Units','normalized', ...
    'Position',[x0 0.010 w max(top-0.015,0.05)], ...
    'String','Load a data file to start.', ...
    'HorizontalAlignment','left','BackgroundColor',bg, ...
    'ForegroundColor',[0.15 0.15 0.15]);

% ---------------- axes, sliders and the divider ----------------
% Positions are all set by layoutTab, which shares the vertical space
% between the two plots according to s.split.
s.ax1 = axes('Parent',tb,'Units','normalized','Position',[0.330 0.700 0.640 0.265]);
s.ax2 = axes('Parent',tb,'Units','normalized','Position',[0.330 0.080 0.640 0.370]);

s.divAx = axes('Parent',tb,'Units','normalized','Position',[0.330 0.480 0.640 0.012], ...
    'XTick',[],'YTick',[],'XLim',[0 1],'YLim',[0 1],'Box','on', ...
    'Color',[0.86 0.86 0.86],'XColor',[0.6 0.6 0.6],'YColor',[0.6 0.6 0.6], ...
    'ButtonDownFcn',@(a,b)startSplit(tb));
line([0.49 0.50 0.51],[0.5 0.5 0.5],'Parent',s.divAx,'LineStyle','none', ...
    'Marker','.','MarkerSize',6,'Color',[0.35 0.35 0.35], ...
    'HitTest','off','PickableParts','none');       % a grip to aim at

snam = {'pre 1','pre 2','post 1','post 2'};
scol = {GRN,GRN,BLU,BLU};
s.sld    = gobjects(1,4);
s.sldLbl = gobjects(1,4);
for k = 1:4
    s.sldLbl(k) = uicontrol('Parent',tb,'Style','text','Units','normalized', ...
        'Position',[0.250 0.5 0.072 0.024],'String',snam{k}, ...
        'HorizontalAlignment','right','BackgroundColor',bg, ...
        'ForegroundColor',scol{k},'FontWeight','bold');
    % plain 0..1 position; Min/Max are set once and never touched again,
    % the energy is worked out from the scan limits when the knob moves
    s.sld(k) = uicontrol('Parent',tb,'Style','slider','Units','normalized', ...
        'Position',[0.330 0.5 0.640 0.030], ...
        'Min',0,'Max',1,'Value',(k-1)/3,'SliderStep',[0.002 0.02], ...
        'Interruptible','off','BusyAction','cancel','Enable','off', ...
        'Callback',@(a,b)onSliderMoved(tb,k));
    installLiveListener(tb,s.sld(k),k);
end

xlabel(s.ax1,'Energy, [eV]');  ylabel(s.ax1,'raw signal');
title(s.ax1,'Raw data and polynomial background');
xlabel(s.ax2,'Energy, [eV]');  ylabel(s.ax2,'\mu, [cm^2/g]');
title(s.ax2,'Normalized absorption');
grid(s.ax1,'on');  grid(s.ax2,'on');

% ---------------- per-graph zoom and autoscale ----------------
% One set for each graph, so they are completely independent.
s.btnZoom = gobjects(2,3);
zlab = {'+','-','A'};
ztip = {'zoom in','zoom out','autoscale to the data'};
for k = 1:2
    for j = 1:3
        if j == 3
            cb = @(a,b)autoScale(tb,k);
        else
            cb = @(a,b)zoomAxes(tb,k,1.45^(2*j-3));   % j=1 in, j=2 out
        end
        s.btnZoom(k,j) = uicontrol('Parent',tb,'Style','pushbutton', ...
            'Units','normalized','Position',[0.253 0.5 0.068 0.028], ...
            'String',zlab{j},'TooltipString',ztip{j},'Callback',cb);
    end
end

% ---------------- hover readout ----------------
s.txHover = uicontrol('Parent',tb,'Style','text','Units','normalized', ...
    'Position',[0.330 0.966 0.455 0.026],'String','', ...
    'HorizontalAlignment','left','BackgroundColor',bg, ...
    'ForegroundColor',[0.25 0.25 0.25],'FontName','FixedWidth');

% marker and label that follow the pointer.  Made once, never deleted:
% they carry a different tag from the plotted data, so a redraw leaves
% them alone, and NaN data keeps them out of the way until first used.
s.hMark = gobjects(1,2);  s.hTxt = gobjects(1,2);
axl = [s.ax1 s.ax2];
for k = 1:2
    s.hMark(k) = line('Parent',axl(k),'XData',NaN,'YData',NaN, ...
        'LineStyle','none','Marker','o','MarkerSize',8,'LineWidth',1.5, ...
        'MarkerEdgeColor',[0.85 0.15 0.15],'Tag','n2phover', ...
        'HandleVisibility','off','HitTest','off','PickableParts','none');
    s.hTxt(k) = text('Parent',axl(k),'Position',[0 0 0],'String','', ...
        'FontSize',8,'Margin',2,'Clipping','on','Visible','off', ...
        'BackgroundColor',[1 1 0.85],'EdgeColor',[0.7 0.7 0.7], ...
        'Tag','n2phover','HandleVisibility','off','HitTest','off', ...
        'PickableParts','none');
end

set(s.ax1,'ButtonDownFcn',@(a,b)startBox(tb,1));
set(s.ax2,'ButtonDownFcn',@(a,b)startBox(tb,2));

s.split = 0.5;              % half the plot area to each graph
layoutTab(s);
set(pnl,'Visible','on');    % one redraw instead of fifty
end

% ======================================================================
%                    PLOT LAYOUT AND THE DIVIDER
% ======================================================================
function [L,W,regTop,regBot] = plotRegion()
% the rectangle the two graphs share, in normalized tab coordinates
L = 0.330;  W = 0.640;  regTop = 0.965;  regBot = 0.080;
end

% ----------------------------------------------------------------------
function layoutTab(s)
%LAYOUTTAB  Share the vertical space between the two graphs.
%   s.split is the fraction going to the top block (graph + sliders).
[L,W,regTop,regBot] = plotRegion();
H  = regTop - regBot;
sp = min(max(s.split,0.35),0.80);

pitch   = 0.036;  sh = 0.030;      % slider pitch and height
sliderH = 4*pitch;
topBot  = regTop - sp*H;           % bottom of the top block

ax1h = max(sp*H - sliderH - 0.055, 0.06);
set(s.ax1,'Position',[L regTop-ax1h W ax1h]);

for k = 1:4
    yk = topBot + (4-k)*pitch + 0.006;
    set(s.sld(k),   'Position',[L      yk       W     sh]);
    set(s.sldLbl(k),'Position',[0.250  yk+0.001 0.072 0.024]);
end

set(s.divAx,'Position',[L topBot-0.020 W 0.012]);

ax2h = max(topBot - 0.020 - 0.045 - regBot, 0.06);
set(s.ax2,'Position',[L regBot W ax2h]);

% Zoom buttons sit in the margin ABOVE each graph, at the right hand
% end: to the left of the axes is where the y tick labels and the y
% label are drawn, so nothing may go there.
axTop = [regTop+0.001, regBot+ax2h+0.008];
for k = 1:2
    for j = 1:3
        set(s.btnZoom(k,j), ...
            'Position',[0.795+(j-1)*0.060 axTop(k) 0.055 0.026]);
    end
end
end

% ----------------------------------------------------------------------
function startSplit(tb)
s = st(tb);
set(s.fig,'WindowButtonMotionFcn',@(a,b)doSplit(tb), ...
          'WindowButtonUpFcn',@(a,b)endSplit(tb),'Pointer','crosshair');
end

function doSplit(tb)
s = st(tb);
if ~ishghandle(s.divAx), return; end
% CurrentPoint is in the divider's own data units, which span its height,
% so this stays correct whatever the tab and figure are doing
cp = get(s.divAx,'CurrentPoint');
p  = get(s.divAx,'Position');
[~,~,regTop,regBot] = plotRegion();
dy = (cp(1,2)-0.5)*p(4);
s.split = min(max(s.split - dy/(regTop-regBot),0.35),0.80);
stSet(tb,s);
layoutTab(s);
end

function endSplit(tb)
s = st(tb);
armHover(s.fig);
end

% ----------------------------------------------------------------------
function resetSplit(fig)
tg = getappdata(fig,'tabgroup');
tb = get(tg,'SelectedTab');
if isempty(tb), return; end
s = st(tb);  s.split = 0.5;  stSet(tb,s);
layoutTab(s);
end

% ----------------------------------------------------------------------
function popOut(fig,which)
%POPOUT  Copy one graph into an ordinary figure window, which can then be
%   resized, zoomed, printed or saved with the usual MATLAB tools.
tg = getappdata(fig,'tabgroup');
tb = get(tg,'SelectedTab');
if isempty(tb), return; end
s = st(tb);
if which == 1, src = s.ax1; else src = s.ax2; end

f  = figure('Name',[get(tb,'Title') ' - ' get(get(src,'Title'),'String')], ...
    'NumberTitle','off','Color','w');   % toolbar left on: zoom, print, save
ax = copyobj(src,f);
set(ax,'Units','normalized','Position',[0.13 0.12 0.80 0.80]);
% the copies must not drive the original tab
set(ax,'ButtonDownFcn','');
set(findall(ax,'Type','line'),'ButtonDownFcn','','HitTest','off');
try, legend(ax,'Location','best'); catch, end
end

% ----------------------------------------------------------------------
function head(pnl,pos,str,bg)
uicontrol(pnl,'Style','text','Units','normalized','Position',pos, ...
    'String',str,'HorizontalAlignment','left','FontWeight','bold', ...
    'BackgroundColor',bg);
end

function lbl(pnl,pos,str,bg)
uicontrol(pnl,'Style','text','Units','normalized','Position',pos, ...
    'String',str,'HorizontalAlignment','left','BackgroundColor',bg);
end

function ed = pairEdit(pnl,pos,str,bg,col,cb)
uicontrol(pnl,'Style','text','Units','normalized', ...
    'Position',[pos(1) pos(2)+0.004 0.55*pos(3) pos(4)-0.004],'String',str, ...
    'HorizontalAlignment','left','BackgroundColor',bg,'ForegroundColor',col);
ed = uicontrol(pnl,'Style','edit','Units','normalized', ...
    'Position',[pos(1)+0.57*pos(3) pos(2) 0.43*pos(3) pos(4)], ...
    'String','','BackgroundColor','w','ForegroundColor',col,'Callback',cb);
end

% ======================================================================
%                             CALLBACKS
% ======================================================================
function loadFile(tb,fname)
s = st(tb);
if nargin < 2
    [f,p] = uigetfile({ ...
        '*.txt;*.dat;*.tbl;*.TBL;*.xmu;*.chi;*.nor;*.mat','Data files'; ...
        '*.*','All files'},'Select a scan file',s.path);
    if isequal(f,0), return; end
else
    [p,b,e] = fileparts(fname);
    if isempty(p), p = pwd; end
    f = [b e];
end

try
    M = readAnyTable(fullfile(p,f));
catch err
    errordlg(sprintf('Could not read %s:\n%s',f,err.message),'Read error');
    return
end
if size(M,2) < 2
    errordlg('The file needs at least two numeric columns.','Read error');
    return
end

s.data = M;  s.file = f;  s.path = p;
s.result = [];  s.stale = true;
set(s.txFile,'String',sprintf('%s   (%d x %d)',f,size(M,1),size(M,2)));
cols = cell(1,size(M,2));
for k = 1:size(M,2), cols{k} = sprintf('%d',k); end
% Value first: a popupmenu may never hold a Value past the end of its
% String, so loading a narrower file after a wide one would otherwise
% error before the new String is applied.
nc = size(M,2);
ec = min(max(1,s.wantCols(1)),nc);          % a duplicate keeps its columns
yc = min(max(1,s.wantCols(2)),nc);
if yc == ec, yc = min(ec+1,nc); end
set(s.popEcol,'Value',1);  set(s.popEcol,'String',cols);
set(s.popYcol,'Value',1);  set(s.popYcol,'String',cols);
set(s.popEcol,'Value',ec);
set(s.popYcol,'Value',yc);
stSet(tb,s);

set(tb,'Title',tabName(f));                 % name the tab after the file
onColumns(tb,true);
end

% ----------------------------------------------------------------------
function onColumns(tb,resetRange)
if nargin < 2, resetRange = false; end
s = st(tb);
if isempty(s.data), return; end

x = s.data(:,get(s.popEcol,'Value'));
y = s.data(:,get(s.popYcol,'Value'));
ok = isfinite(x) & isfinite(y);
x = x(ok);  y = y(ok);
[x,ix] = sort(x);  y = y(ix);

switch get(s.popUnit,'Value')
    case 1                                  % eV, the default: as read
    case 2                                  % keV
        x = x*1000;
    case 3                                  % auto
        if max(x) < 1000, x = x*1000; end
end

s.x = x;  s.y = y;
s.result = [];  s.stale = true;
s.zoomed = [false false];        % a new scan starts from the full view

% a duplicated tab keeps the windows it was given, once, if they fit
r = getRange(s);
carry = s.keepRange && all(isfinite(r)) && all(diff(r) > 0) ...
        && r(2) > x(1) && r(3) < x(end);
if carry
    r(1) = max(r(1),x(1));  r(4) = min(r(4),x(end));
    s.keepRange = false;
end
stSet(tb,s);

if carry
    setRange(s,r);
    afterParamChange(tb);
    return
end

needRange = resetRange || any(~isfinite(getRange(s)));
if needRange && ~isempty(s.Z)
    defaultRange(tb);                       % also replots
else
    refreshPlots(tb);
    s = st(tb);
    syncSliders(s);
    if isempty(s.Z)
        setStatus(s,{sprintf('%d points, %.1f to %.1f eV',numel(x),x(1),x(end)), ...
                     'Now give the element and edge.'});
    else
        afterParamChange(tb);
    end
end
end

% ----------------------------------------------------------------------
function onElement(tb)
s   = st(tb);
el  = elements(s);
str = strtrim(get(s.edZ,'String'));
if isempty(str), return; end

Z = str2double(str);
if isnan(Z)
    try
        Z = zymbol(str);
    catch
        Z = [];
    end
end
if isempty(Z) || numel(Z)~=1 || ~isfinite(Z) || Z<1 || Z>94 || Z~=fix(Z)
    errordlg(sprintf('"%s" is not a usable element (Z = 1 ... 94).',str),'Element');
    return
end
if el(Z).K_edge <= 0
    errordlg(sprintf('McMaster has no data for %s (Z = %d).',el(Z).symbol,Z),'Element');
    return
end

s.Z = Z;
set(s.edZ,'String',el(Z).symbol);
stSet(tb,s);
onEdge(tb);
end

% ----------------------------------------------------------------------
function onEdge(tb)
s  = st(tb);
el = elements(s);
if isempty(s.Z), return; end

k = get(s.popEdge,'Value');
if k == 1
    E0 = el(s.Z).K_edge*1000;   name = 'K';
else
    E0 = el(s.Z).L_edge(k-1)*1000;   name = sprintf('L%d',k-1);
end
if E0 <= 0
    errordlg(sprintf('McMaster has no %s edge for %s.',name,el(s.Z).symbol),'Edge');
    return
end

set(s.edE0,'String',sprintf('%.1f',E0));
if ~isempty(s.x) && any(~isfinite(getRange(s)))
    defaultRange(tb);
else
    afterParamChange(tb);
end
end

% ----------------------------------------------------------------------
function onE0(tb)
s = st(tb);
if ~isfinite(getE0(s))
    setStatus(s,{'Edge energy must be a number, in eV.'});
    return
end
afterParamChange(tb);
end

% ----------------------------------------------------------------------
function defaultRange(tb)
% same defaults as norm2poly.m, around the selected edge:
% scan start .. E0-20 eV,  E0+80 eV .. scan end
s  = st(tb);
E0 = getE0(s);
if isempty(s.x) || ~isfinite(E0)
    setStatus(s,{'Load a file and set the edge energy first.'});
    return
end
r = [s.x(1), E0-20, E0+80, s.x(end)];
setRange(s,r);
afterParamChange(tb);
if any(diff(r) <= 0)
    setStatus(s,{'The edge energy sits outside the scan;', ...
                 'set the four window energies by hand.'});
end
end

% ----------------------------------------------------------------------
function onRange(tb)
s = st(tb);
r = getRange(s);
if any(~isfinite(r))
    setStatus(s,{'Fit windows: please enter four numbers.'});  return
end
syncSliders(s);
afterParamChange(tb);
end

% ----------------------------------------------------------------------
function onAuto(tb)
s = st(tb);
if get(s.chkAuto,'Value'), runFit(tb,false); end
end

% ----------------------------------------------------------------------
function afterParamChange(tb)
% something that affects the fit changed: redraw the raw view, then
% either re-run straight away or wait for the Run button
s = st(tb);
s.stale = true;
stSet(tb,s);
refreshPlots(tb);

s = st(tb);
if get(s.chkAuto,'Value')
    runFit(tb,false);
elseif ~isempty(s.x)
    r = getRange(s);
    if any(diff(r) <= 0)
        setStatus(s,{'The four energies must increase:', ...
            sprintf('%.1f  %.1f  %.1f  %.1f',r(1),r(2),r(3),r(4))});
    else
        setStatus(s,{'Parameters changed.','Press "Run normalization".'});
    end
end
end

% ======================================================================
%                            RUN THE FIT
% ======================================================================
function runFit(tb,loud)
if nargin < 2, loud = false; end
s = st(tb);

if isempty(s.x)
    if loud, errordlg('Load a data file first.','Normalize'); end
    return
end
if isempty(s.Z)
    if loud, errordlg('Give the element first.','Normalize'); end
    return
end
r = getRange(s);
if any(~isfinite(r)) || any(diff(r) <= 0)
    if loud, errordlg('The four fit-window energies must increase.','Normalize'); end
    return
end

n_th = get(s.popOrder,'Value') - 1;
try
    res = n2p_compute(s.x,s.y,s.Z,n_th,r,elements(s));
catch err
    s.result = [];  s.stale = true;  stSet(tb,s);
    refreshPlots(tb);
    setStatus(s,{'Fit failed:',err.message});
    if loud, errordlg(err.message,'Fit failed'); end
    return
end

s.result = res;  s.stale = false;
stSet(tb,s);
refreshPlots(tb);

s = st(tb);
el = elements(s);
setStatus(s,{ ...
    sprintf('%s, order %d',el(s.Z).symbol,n_th), ...
    sprintf('%d pre + %d post points',res.n1,res.n2), ...
    sprintf('McMaster step: %.1f cm^2/g',res.mu_ref(end)), ...
    sprintf('fit residual rms: %.3g',res.rms)});
end

% ======================================================================
%                            COMPUTATION
% ======================================================================
function res = n2p_compute(x,y,Z,n_th,r,element)
%N2P_COMPUTE  the maths of norm2poly.m, with no graphics attached.
%   x    energy [eV], ascending        r     [pre1 pre2 post1 post2]
%   y    raw signal                    n_th  polynomial order

idx_pre  = x >= r(1) & x <= r(2);
idx_post = x >= r(3) & x <= r(4);
n1 = sum(idx_pre);  n2 = sum(idx_post);

if n1 < 2 || n2 < 2
    error('each window needs 2+ points (got %d and %d)',n1,n2);
end
if n1+n2 < n_th+3
    error('%d points is too few for order %d',n1+n2,n_th);
end
if x(1)/1000 < element(Z).M_edge
    error('scan starts below the M edge (%.1f eV)',element(Z).M_edge*1000);
end

x_1 = x(idx_pre);   x_2 = x(idx_post);   x_12 = [x_1; x_2];
y_1 = y(idx_pre);   y_2 = y(idx_post);   y_12 = [y_1; y_2];

% --- McMaster cross-section, pre-edge baseline extrapolated away ------
mu = McMaster(Z,x_12/1000);                  % [cm^2/g], energies in keV
p  = polyfit(x_1,mu(1:n1),1);
mu = mu - [mu(1:n1); polyval(p,x_2)];        % identically zero below the edge

% --- map onto [-1,1] -------------------------------------------------
[Nx_12,mid_x,span_x]   = unify(x_12);
[Ny_12,mid_y,span_y]   = unify(y_12);
[Nmu,  mid_mu,span_mu] = unify(mu);

% --- design matrix: edge step + Legendre polynomials ------------------
A = [Nmu, ones(numel(x_12),1)];
if n_th >= 1, A = [A, Nx_12]; end
for j = 2:n_th
    A = [A, polyval(legpoly(j),Nx_12)]; %#ok<AGROW>
end
C = A\Ny_12(:);
if C(1) == 0, error('the edge-step coefficient came out zero'); end

% --- apply over the whole scan ---------------------------------------
Nx = (x-mid_x)/span_x;
Ny = (y-mid_y)/span_y;
T  = ones(numel(x),1);
if n_th >= 1, T = [T, Nx]; end
for j = 2:n_th
    T = [T, polyval(legpoly(j),Nx)]; %#ok<AGROW>
end
bg   = T*C(2:end);                            % background, scaled units
Xabs = (Ny-bg)*span_mu/C(1) + mid_mu;         % normalized signal [cm^2/g]

% --- McMaster reference over the whole scan ---------------------------
mu_ref = McMaster(Z,x/1000);
below  = x < r(2);
if sum(below) >= 2
    q = polyfit(x(below),mu_ref(below),1);
    mu_ref = mu_ref - polyval(q,x);
    mu_ref(below) = 0;
else
    mu_ref = mu_ref - mu_ref(1);
end

res = struct( ...
    'x',x, 'y',y, 'Xabs',Xabs, 'mu_ref',mu_ref, ...
    'bg_raw',bg*span_y+mid_y, ...             % background in raw units
    'C',C, 'n_th',n_th, 'range',r, 'n1',n1, 'n2',n2, ...
    'idx_pre',idx_pre, 'idx_post',idx_post, ...
    'rms',sqrt(mean((y_12 - (A*C)*span_y - mid_y).^2)));
end

% ======================================================================
%                             PLOTTING
% ======================================================================
function clearPlot(ax)
% Delete everything this GUI drew.  cla() leaves objects whose
% HandleVisibility is 'off' behind - which is why the dashed limits used
% to pile up - so they are tagged and deleted explicitly instead.
delete(findall(ax,'Tag','n2p'));
end

% ----------------------------------------------------------------------
function refreshPlots(tb)
s = st(tb);
if isempty(s.x)
    clearPlot(s.ax1);  clearPlot(s.ax2);  return
end
s = drawAx1(s);
s = drawAx2(s);
stSet(tb,s);
end

% ----------------------------------------------------------------------
function s = drawAx1(s)
% raw scan, shaded fit windows, draggable limits and, once the fit has
% been run, the fitted background
[xl,yl] = viewLimits(s,1);      % keeps a hand-set zoom across redraws
r  = getRange(s);

clearPlot(s.ax1);
hold(s.ax1,'on');

s.patch = gobjects(1,2);
j   = [1 3];
col = {s.colPre, s.colPost};
for k = 1:2
    if all(isfinite(r([j(k) j(k)+1]))) && r(j(k)+1) > r(j(k))
        s.patch(k) = patch('Parent',s.ax1,'Tag','n2p', ...
            'XData',[r(j(k)) r(j(k)+1) r(j(k)+1) r(j(k))], ...
            'YData',[yl(1) yl(1) yl(2) yl(2)], ...
            'FaceColor',col{k},'FaceAlpha',0.12,'EdgeColor','none', ...
            'HitTest','off','PickableParts','none','HandleVisibility','off');
    end
end

hRaw = plot(s.ax1,s.x,s.y,'k.-','MarkerSize',6,'Tag','n2p', ...
    'DisplayName','raw data','HitTest','off','PickableParts','none');
s.curveH{1} = hRaw;  s.curveName{1} = {'raw'};
if ~isempty(s.result) && ~s.stale
    hBg = plot(s.ax1,s.result.x,s.result.bg_raw,'r-','LineWidth',1.5, ...
        'Tag','n2p','DisplayName','background', ...
        'HitTest','off','PickableParts','none');
    s.curveH{1}(end+1) = hBg;  s.curveName{1}{end+1} = 'background';
end

set(s.ax1,'XLim',xl,'YLim',yl);
xlabel(s.ax1,'Energy, [eV]');  ylabel(s.ax1,'raw signal');
title(s.ax1,sprintf('%s   -   raw data',strrep(s.file,'_','\_')));
grid(s.ax1,'on');
legend(s.ax1,'Location','best');

s.line = gobjects(1,4);
lcol = {s.colPre,s.colPre,s.colPost,s.colPost};
for k = 1:4
    if isfinite(r(k))
        s.line(k) = plot(s.ax1,r(k)*[1 1],yl,'--','Color',lcol{k}, ...
            'LineWidth',2,'Tag','n2p','HandleVisibility','off', ...
            'ButtonDownFcn',@(a,b)startDrag(s.tab,k));
    end
end
hold(s.ax1,'off');
end

% ----------------------------------------------------------------------
function s = drawAx2(s)
clearPlot(s.ax2);
s.curveH{2} = gobjects(1,0);  s.curveName{2} = {};
xlabel(s.ax2,'Energy, [eV]');
ylabel(s.ax2,'X-ray absorption coefficient, [cm^2/g]');
grid(s.ax2,'on');
title(s.ax2,'Normalized absorption');

if isempty(s.result) || s.stale
    [xl0,yl0] = viewLimits(s,2);
    set(s.ax2,'XLim',xl0,'YLim',yl0);
    text(0.5,0.5,'press  "Run normalization"','Parent',s.ax2, ...
        'Units','normalized','HorizontalAlignment','center','Tag','n2p', ...
        'Color',[0.5 0.5 0.5],'FontAngle','italic');
    return
end

r = s.result;
hold(s.ax2,'on');
plot(s.ax2,[r.x(1) r.x(end)],[0 0],'k:','Tag','n2p', ...
    'HandleVisibility','off','HitTest','off','PickableParts','none');
h1 = plot(s.ax2,r.x,r.Xabs,'b-','LineWidth',1.2,'Tag','n2p', ...
    'DisplayName','normalized','HitTest','off','PickableParts','none');
h2 = plot(s.ax2,r.x,r.mu_ref,'c-','LineWidth',1.5,'Tag','n2p', ...
    'DisplayName','McMaster','HitTest','off','PickableParts','none');
s.curveH{2} = [h1 h2];  s.curveName{2} = {'normalized','McMaster'};
[xl2,yl2] = viewLimits(s,2);
set(s.ax2,'XLim',xl2,'YLim',yl2);
legend(s.ax2,'Location','southeast');
hold(s.ax2,'off');
end

% ----------------------------------------------------------------------
function v = padded(v)
if ~all(isfinite(v)), v = [0 1]; return, end
if v(2) <= v(1)
    v = v(1) + [-1 1]*max(abs(v(1))*0.05,1);
else
    v = v + 0.05*(v(2)-v(1))*[-1 1];
end
end

% ======================================================================
%                    ZOOM, AUTOSCALE AND VIEW LIMITS
% ======================================================================
function [xl,yl] = autoLimits(s,k)
%AUTOLIMITS  The natural full view of one graph.
if k == 1
    xl = [s.x(1) s.x(end)];
    yl = padded([min(s.y) max(s.y)]);
elseif isempty(s.result) || s.stale
    xl = [s.x(1) s.x(end)];
    yl = [0 1];
else
    r  = s.result;
    xl = [r.x(1) r.x(end)];
    yl = padded([min([r.Xabs;r.mu_ref]) max([r.Xabs;r.mu_ref])]);
end
if xl(2) <= xl(1), xl = xl(1) + [0 1]; end
end

% ----------------------------------------------------------------------
function [xl,yl] = viewLimits(s,k)
%VIEWLIMITS  What a redraw should use: the hand-set zoom if there is one,
%   otherwise the full view.  This is what stops the zoom being thrown
%   away every time the fit is recomputed.
ax = axesOf(s,k);
if s.zoomed(k) && ishghandle(ax)
    xl = get(ax,'XLim');  yl = get(ax,'YLim');
else
    [xl,yl] = autoLimits(s,k);
end
end

function ax = axesOf(s,k)
if k == 1, ax = s.ax1; else ax = s.ax2; end
end

% ----------------------------------------------------------------------
function zoomAxes(tb,k,f)
%ZOOMAXES  Zoom one graph about the middle of its current view.
s  = st(tb);
ax = axesOf(s,k);
xl = get(ax,'XLim');  yl = get(ax,'YLim');
applyZoom(tb,s,k,mean(xl),mean(yl),f);
end

% ----------------------------------------------------------------------
function applyZoom(tb,s,k,cx,cy,f)
ax = axesOf(s,k);
xl = get(ax,'XLim');  yl = get(ax,'YLim');
nx = cx + (xl-cx)*f;
ny = cy + (yl-cy)*f;
if diff(nx) <= 0 || diff(ny) <= 0, return; end
if diff(nx) > 1e12 || diff(ny) > 1e12, return; end
set(ax,'XLim',nx,'YLim',ny);
s.zoomed(k) = true;
stSet(tb,s);
if k == 1, stretchWindowMarks(s); end
end

% ----------------------------------------------------------------------
function autoScale(tb,k)
%AUTOSCALE  Back to the full view of that graph.
s = st(tb);
s.zoomed(k) = false;
stSet(tb,s);
[xl,yl] = autoLimits(s,k);
set(axesOf(s,k),'XLim',xl,'YLim',yl);
if k == 1, stretchWindowMarks(s); end
end

% ----------------------------------------------------------------------
function stretchWindowMarks(s)
% the shaded windows and the dashed limits are drawn to fill the axes, so
% they have to be restretched whenever the top graph's y range changes
yl = get(s.ax1,'YLim');
for k = 1:4
    if numel(s.line) >= k && ishghandle(s.line(k))
        set(s.line(k),'YData',yl);
    end
end
for k = 1:2
    if numel(s.patch) >= k && ishghandle(s.patch(k))
        set(s.patch(k),'YData',[yl(1) yl(1) yl(2) yl(2)]);
    end
end
end

% ----------------------------------------------------------------------
function onScroll(fig,evt)
%ONSCROLL  Wheel zoom about the pointer, on whichever graph it is over.
tg = getappdata(fig,'tabgroup');
if isempty(tg) || ~ishghandle(tg), return; end
tb = get(tg,'SelectedTab');
if isempty(tb), return; end
s = st(tb);
if isempty(s) || ~isstruct(s) || isempty(s.x), return; end

f = 1.2^double(evt.VerticalScrollCount);     % wheel up = zoom in
for k = 1:2
    [in,cx,cy] = pointerIn(s,k);
    if in, applyZoom(tb,s,k,cx,cy,f); return, end
end
end

% ----------------------------------------------------------------------
function [in,cx,cy] = pointerIn(s,k)
%POINTERIN  Is the pointer inside graph k, and where in data units?
ax = axesOf(s,k);
in = false;  cx = NaN;  cy = NaN;
if ~ishghandle(ax), return; end
cp = get(ax,'CurrentPoint');
xl = get(ax,'XLim');  yl = get(ax,'YLim');
cx = cp(1,1);  cy = cp(1,2);
in = cx >= xl(1) && cx <= xl(2) && cy >= yl(1) && cy <= yl(2);
end

% ======================================================================
%                        BOX ZOOM (RUBBER BAND)
% ======================================================================
function startBox(tb,k)
%STARTBOX  Press anywhere on a graph and drag to outline the region to
%   zoom into.  A press that lands on one of the dashed window limits
%   moves that limit instead - the line gets the click first.
s = st(tb);
if isempty(s) || ~isstruct(s) || isempty(s.x), return; end
if ~strcmp(get(s.fig,'SelectionType'),'normal'), return; end   % left button

[in,cx,cy] = pointerIn(s,k);
if ~in, return; end

ax = axesOf(s,k);
s.boxK = k;
s.boxAnchor = [cx cy];
s.boxH = line('Parent',ax,'XData',cx*ones(1,5),'YData',cy*ones(1,5), ...
    'Color',[0.20 0.20 0.20],'LineStyle','--','LineWidth',1, ...
    'Tag','n2pbox','HandleVisibility','off','HitTest','off', ...
    'PickableParts','none');
stSet(tb,s);

set(s.fig,'WindowButtonMotionFcn',@(a,b)doBox(tb), ...
          'WindowButtonUpFcn',@(a,b)endBox(tb),'Pointer','crosshair');
end

% ----------------------------------------------------------------------
function doBox(tb)
s = st(tb);
if ~s.boxK || isempty(s.boxH) || ~ishghandle(s.boxH), return; end
ax = axesOf(s,s.boxK);
cp = get(ax,'CurrentPoint');
xl = get(ax,'XLim');  yl = get(ax,'YLim');
cx = min(max(cp(1,1),xl(1)),xl(2));       % stay inside the current view
cy = min(max(cp(1,2),yl(1)),yl(2));
a  = s.boxAnchor;
set(s.boxH,'XData',[a(1) cx cx a(1) a(1)],'YData',[a(2) a(2) cy cy a(2)]);
end

% ----------------------------------------------------------------------
function endBox(tb)
s = st(tb);
k = s.boxK;
X = [];  Y = [];
if ~isempty(s.boxH) && ishghandle(s.boxH)
    X = get(s.boxH,'XData');  Y = get(s.boxH,'YData');
    delete(s.boxH);
end
s.boxK = 0;  s.boxH = gobjects(1,0);
stSet(tb,s);
armHover(s.fig);

if k == 0 || isempty(X), return; end
ax = axesOf(s,k);
xl = get(ax,'XLim');  yl = get(ax,'YLim');
nx = [min(X) max(X)];  ny = [min(Y) max(Y)];
% a plain click, or a sliver of a box, is not a zoom request
if diff(nx) < 0.01*diff(xl) || diff(ny) < 0.01*diff(yl), return; end

set(ax,'XLim',nx,'YLim',ny);
s.zoomed(k) = true;
stSet(tb,s);
if k == 1, stretchWindowMarks(s); end
end

% ======================================================================
%                       HOVER READOUT ON THE GRAPHS
% ======================================================================
function armHover(fig)
%ARMHOVER  Mouse motion drives the readout whenever nothing is being
%   dragged.  The drag handlers borrow WindowButtonMotionFcn and hand it
%   back here when they are done.
set(fig,'WindowButtonMotionFcn',@(a,b)onHover(fig), ...
        'WindowButtonUpFcn','','Pointer','arrow');
end

% ----------------------------------------------------------------------
function onHover(fig)
tg = getappdata(fig,'tabgroup');
if isempty(tg) || ~ishghandle(tg), return; end
tb = get(tg,'SelectedTab');
if isempty(tb) || ~ishghandle(tb), return; end
s = st(tb);
if isempty(s) || ~isstruct(s) || isempty(s.x), return; end

for k = 1:2
    [in,cx,cy] = pointerIn(s,k);
    if in, showHover(s,k,cx,cy);  return, end
end
clearHover(s);
end

% ----------------------------------------------------------------------
function showHover(s,k,cx,cy)
%SHOWHOVER  Report the cursor position and the nearest point on whichever
%   curve in that graph the pointer is closest to.  Point at a curve and
%   the readout follows it; move across to another and it switches.

hs = s.curveH{k};  names = s.curveName{k};
if isempty(hs)
    set(s.txHover,'String',sprintf('cursor %.1f, %.4g',cx,cy));
    clearMarks(s);  return
end

% distances measured in fractions of the visible axes, so that x in eV
% and y in cm^2/g can be compared sensibly
ax  = axesOf(s,k);
xl  = get(ax,'XLim');  yl = get(ax,'YLim');
sx  = max(diff(xl),eps);  sy = max(diff(yl),eps);

best = inf;  px = [];  py = [];  ci = 1;  i = 1;
for c = 1:numel(hs)
    if ~ishghandle(hs(c)), continue, end
    X = get(hs(c),'XData');  Y = get(hs(c),'YData');
    if numel(X) < 2, continue, end
    d = ((X(:)-cx)/sx).^2 + ((Y(:)-cy)/sy).^2;
    [dm,im] = min(d);
    if dm < best
        best = dm;  ci = c;  i = im;  px = X(:);  py = Y(:);
    end
end
if isempty(px)
    set(s.txHover,'String',sprintf('cursor %.1f, %.4g',cx,cy));
    clearMarks(s);  return
end

extra = '';
if k == 2 && ~isempty(s.result) && ~s.stale && numel(s.result.mu_ref) == numel(px)
    extra = sprintf('   (McMaster %.5g)',s.result.mu_ref(i));
end

set(s.txHover,'String',sprintf('cursor %.1f, %.4g   |   %s %.1f, %.5g%s', ...
    cx,cy,names{ci},px(i),py(i),extra));

other = 3-k;
set(s.hMark(other),'XData',NaN,'YData',NaN);
set(s.hTxt(other),'Visible','off');

set(s.hMark(k),'XData',px(i),'YData',py(i), ...
    'MarkerEdgeColor',get(hs(ci),'Color'),'MarkerFaceColor','none');

% keep the label inside the axes: flip it to the left near the right edge
xl = get(get(s.hMark(k),'Parent'),'XLim');
if px(i) > mean(xl)
    ha = 'right';  dx = -0.01*(xl(2)-xl(1));
else
    ha = 'left';   dx =  0.01*(xl(2)-xl(1));
end
set(s.hTxt(k),'Position',[px(i)+dx py(i) 0], ...
    'String',sprintf(' %s:  %.1f eV,  %.5g ',names{ci},px(i),py(i)), ...
    'HorizontalAlignment',ha,'VerticalAlignment','bottom','Visible','on');
end

% ----------------------------------------------------------------------
function clearMarks(s)
for k = 1:2
    if ishghandle(s.hMark(k)), set(s.hMark(k),'XData',NaN,'YData',NaN); end
    if ishghandle(s.hTxt(k)),  set(s.hTxt(k),'Visible','off'); end
end
end

function clearHover(s)
set(s.txHover,'String','');
clearMarks(s);
end

% ======================================================================
%                        THE FOUR WINDOW SLIDERS
% ======================================================================
function installLiveListener(tb,hs,k)
% live updates while the knob moves; fall back to release-only if this
% MATLAB has neither event
try
    addlistener(hs,'ContinuousValueChange',@(a,b)onSlider(tb,k));
    return
catch
end
try
    addlistener(hs,'Value','PostSet',@(a,b)onSlider(tb,k));
catch
end
end

% ----------------------------------------------------------------------
function onSlider(tb,k)
% Writes ONLY its own edit box and moves the existing dashed line; it
% never redraws the axes and never writes back to a slider.
if getappdata(tb,'lock'), return; end
s = st(tb);
if isempty(s) || ~isstruct(s) || isempty(s.x), return; end   % still building

lo = s.x(1);  hi = s.x(end);
E  = lo + min(max(get(s.sld(k),'Value'),0),1)*(hi-lo);

eds = [s.edPre1 s.edPre2 s.edPost1 s.edPost2];
set(eds(k),'String',sprintf('%.1f',E));

updateWindowGraphics(s);
r = getRange(s);
setStatus(s,{sprintf('pre-edge    %.1f  to  %.1f',r(1),r(2)), ...
             sprintf('post-edge  %.1f  to  %.1f',r(3),r(4)),''});
end

% ----------------------------------------------------------------------
function onSliderMoved(tb,k)
onSlider(tb,k);
afterParamChange(tb);
end

% ======================================================================
%                    DRAGGING THE WINDOW BOUNDARIES
% ======================================================================
function startDrag(tb,k)
s = st(tb);
s.drag = k;
stSet(tb,s);
set(s.fig,'WindowButtonMotionFcn',@(a,b)doDrag(tb), ...
          'WindowButtonUpFcn',@(a,b)endDrag(tb),'Pointer','fleur');
end

function doDrag(tb)
s = st(tb);
if ~s.drag, return; end
cp = get(s.ax1,'CurrentPoint');
xv = cp(1,1);
if ~isfinite(xv), return; end
r  = getRange(s);
% held inside the scan, but free to pass its neighbours; the ordering is
% only checked when the fit is run
r(s.drag) = min(max(xv,s.x(1)),s.x(end));
setRange(s,r);
end

function endDrag(tb)
s = st(tb);
s.drag = 0;
stSet(tb,s);
armHover(s.fig);            % motion goes back to the readout
afterParamChange(tb);
end

% ======================================================================
%                              OUTPUT
% ======================================================================
function saveNrm(tb)
s = st(tb);
if isempty(s.result) || s.stale
    errordlg('Run the normalization first.','Save');  return
end
base = strtok(s.file,'.');
if isempty(base), base = 'output'; end
[f,p] = uiputfile({'*.nrm','normalized data (*.nrm)';'*.*','All files'}, ...
    'Save normalized data',fullfile(s.path,[base '.nrm']));
if isequal(f,0), return; end

out = [s.result.x, s.result.Xabs, s.result.mu_ref]; %#ok<NASGU>
try
    save(fullfile(p,f),'out','-ascii');
    setStatus(s,{'saved:',f,'columns: energy, normalized, McMaster'});
catch err
    errordlg(err.message,'Save failed');
end
end

% ----------------------------------------------------------------------
function toWorkspace(tb)
s = st(tb);
if isempty(s.result) || s.stale
    errordlg('Run the normalization first.','Export');  return
end
name = varName(s.file);
assignin('base',name,s.result);
setStatus(s,{['struct "' name '" is now'],'in the base workspace'});
end

% ======================================================================
%                             HELPERS
% ======================================================================
function r = getRange(s)
r = [str2double(get(s.edPre1,'String')),  str2double(get(s.edPre2,'String')), ...
     str2double(get(s.edPost1,'String')), str2double(get(s.edPost2,'String'))];
end

function setRange(s,r)
set(s.edPre1, 'String',sprintf('%.1f',r(1)));
set(s.edPre2, 'String',sprintf('%.1f',r(2)));
set(s.edPost1,'String',sprintf('%.1f',r(3)));
set(s.edPost2,'String',sprintf('%.1f',r(4)));
syncSliders(s);
updateWindowGraphics(s);
end

function E0 = getE0(s)
E0 = str2double(get(s.edE0,'String'));
end

function setStatus(s,lines)
set(s.txStat,'String',lines);
end

% ----------------------------------------------------------------------
function syncSliders(s)
% knobs follow the edit boxes.  Min/Max are never touched, and the lock
% keeps these programmatic writes out of onSlider.
if ~isfield(s,'sld'), return; end
if isempty(s.x)
    set(s.sld,'Enable','off');  return
end

lo = s.x(1);  hi = s.x(end);
if hi <= lo, hi = lo + 1; end
r = getRange(s);

setappdata(s.tab,'lock',1);
for k = 1:4
    p = (r(k)-lo)/(hi-lo);
    if ~isfinite(p), p = (k-1)/3; end
    set(s.sld(k),'Value',min(max(p,0),1));
end
setappdata(s.tab,'lock',0);

set(s.sld,'Enable','on');
end

% ----------------------------------------------------------------------
function updateWindowGraphics(s)
% move the existing dashed limits and shaded windows in place - no new
% objects are ever created here
r = getRange(s);

for k = 1:4
    if isfield(s,'line') && numel(s.line) >= k && ishghandle(s.line(k)) ...
            && isfinite(r(k))
        set(s.line(k),'XData',[r(k) r(k)]);
    end
end

j = [1 3];
for k = 1:2
    if isfield(s,'patch') && numel(s.patch) >= k && ishghandle(s.patch(k)) ...
            && all(isfinite(r([j(k) j(k)+1])))
        set(s.patch(k),'XData',[r(j(k)) r(j(k)+1) r(j(k)+1) r(j(k))]);
    end
end
end

% ----------------------------------------------------------------------
function t = tabName(f)
t = strtok(f,'.');
if isempty(t), t = 'dataset'; end
if numel(t) > 16, t = [t(1:15) '~']; end
end

function v = varName(f)
v = regexprep(strtok(f,'.'),'\W','_');
if isempty(v) || ~isletter(v(1)), v = ['n2p_' v]; end
if numel(v) > 40, v = v(1:40); end
end

% ----------------------------------------------------------------------
function M = readAnyTable(full)
%READANYTABLE  numeric table from an ascii file (headers tolerated) or .mat
[~,~,ext] = fileparts(full);
if strcmpi(ext,'.mat')
    S  = load(full);
    fn = fieldnames(S);
    M  = [];
    for k = 1:numel(fn)
        v = S.(fn{k});
        if isnumeric(v) && size(v,2) >= 2 && size(v,1) > 2
            M = v;  break
        end
    end
    if isempty(M)
        error('no numeric two-column variable found in the file');
    end
    return
end

d = importdata(full);
if isstruct(d)
    M = d.data;
elseif isnumeric(d)
    M = d;
else
    error('unrecognized file contents');
end
if isempty(M) || ~isnumeric(M)
    error('no numeric data found');
end
M = M(all(isfinite(M(:,1:min(2,size(M,2)))),2),:);
end

% ======================================================================
%      EMBEDDED McMASTER TABLES   (the former element.mat, in code)
% ======================================================================
function el = n2p_elementTable()
%N2P_ELEMENTTABLE  The contents of element.mat, written out as source so
%   that this GUI is a single self-contained file.  Values are byte-for-
%   byte the ones in element.mat (May 1969 McMaster compilation, as
%   rebuilt by make_element.m in 2005).  Built once, then cached.
%
%   Po, At, Fr, Ra, Ac, Pa and Np are all zeros: McMaster did not
%   tabulate them.

persistent EL
if ~isempty(EL), el = EL;  return, end

sym = { ...
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', ...
    'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', ...
    'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', ...
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', ...
    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', ...
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', ...
    'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', ...
    'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu' };

% scalars, one row per element:
%   weight  density  m2a  K_edge  M_edge  K_alpha1  K_beta1  L_alpha1  L_beta1  fluo_K
A = [ ...
    1.008, 8.9e-05, 1.674, 0.014, 0, 0, 0, 0, 0, 0
    4.003, 0.000179, 6.647, 0.025, 0, 0, 0, 0, 0, 0
    6.94, 0.534, 11.52, 0.055, 0, 0, 0, 0, 0, 0
    9.012, 1.848, 14.96, 0.112, 0, 0, 0, 0, 0, 0
    10.811, 2.34, 17.95, 0.188, 0, 0, 0, 0, 0, 0.0017
    12.01, 2.26, 19.94, 0.295, 0, 0, 0, 0, 0, 0.0028
    14.008, 0.00125, 23.26, 0.402, 0, 0, 0, 0, 0, 0.0052
    16, 0.001429, 26.57, 0.537, 0, 0, 0, 0, 0, 0.0083
    19, 1.108, 31.55, 0.686, 0, 0, 0, 0, 0, 0.013
    20.183, 0.0009, 33.51, 0.867, 0, 0, 0, 0, 0, 0.018
    22.997, 0.97, 38.19, 1.072, 0, 1.041, 1.067, 0, 0, 0.023
    24.32, 1.74, 40.38, 1.305, 0, 1.254, 1.297, 0, 0, 0.03
    26.97, 2.72, 44.78, 1.56, 0, 1.487, 1.553, 0, 0, 0.039
    28.086, 2.33, 46.63, 1.839, 0, 1.74, 1.832, 0, 0, 0.05
    30.975, 1.82, 51.43, 2.149, 0, 2.015, 2.136, 0, 0, 0.063
    32.066, 2, 53.24, 2.472, 0.017, 2.308, 2.464, 0, 0, 0.078
    35.457, 1.56, 58.87, 2.822, 0.017, 2.622, 2.815, 0, 0, 0.097
    39.944, 0.001784, 66.32, 3.202, 0.027, 2.957, 3.192, 0, 0, 0.118
    39.102, 0.862, 64.93, 3.607, 0.034, 3.313, 3.589, 0, 0, 0.14
    40.08, 1.55, 66.55, 4.038, 0.044, 3.691, 4.012, 0, 0, 0.163
    44.96, 2.992, 74.65, 4.493, 0.054, 4.09, 4.46, 0, 0, 0.188
    47.9, 4.54, 79.53, 4.965, 0.059, 4.51, 4.931, 0, 0, 0.214
    50.942, 6.11, 84.59, 5.465, 0.067, 4.952, 5.427, 0, 0, 0.243
    51.996, 7.19, 86.34, 5.989, 0.074, 5.414, 5.946, 0, 0, 0.275
    54.94, 7.42, 91.22, 6.54, 0.084, 5.898, 6.49, 0, 0, 0.308
    55.85, 7.86, 92.74, 7.112, 0.094, 6.403, 7.057, 0, 0, 0.34
    58.933, 8.9, 97.85, 7.709, 0.101, 6.93, 7.649, 0, 0, 0.373
    58.69, 8.9, 97.45, 8.333, 0.113, 7.477, 8.264, 0, 0, 0.406
    63.54, 8.94, 105.5, 8.979, 0.12, 8.047, 8.904, 0, 0, 0.44
    65.38, 7.14, 108.6, 9.659, 0.139, 8.638, 9.571, 1.009, 1.032, 0.474
    69.72, 5.903, 115.8, 10.367, 0.158, 9.251, 10.263, 1.096, 1.122, 0.507
    72.59, 5.323, 120.5, 11.104, 0.181, 9.885, 10.981, 1.186, 1.216, 0.535
    74.92, 5.73, 124.4, 11.868, 0.206, 10.543, 11.725, 1.282, 1.317, 0.562
    78.96, 4.79, 131.1, 12.658, 0.232, 11.221, 12.495, 1.419, 1.379, 0.589
    79.92, 3.12, 132.7, 13.474, 0.257, 11.923, 13.29, 1.48, 1.526, 0.618
    83.8, 0.00374, 139.1, 14.322, 0.288, 12.648, 14.112, 1.587, 1.638, 0.643
    85.48, 1.532, 141.9, 15.2, 0.322, 13.394, 14.96, 1.694, 1.752, 0.667
    87.62, 2.54, 145.5, 16.105, 0.358, 14.164, 15.834, 1.806, 1.872, 0.69
    88.905, 4.405, 147.6, 17.08, 0.395, 14.957, 16.736, 1.922, 1.996, 0.71
    91.22, 6.53, 151.5, 17.998, 0.431, 15.774, 17.666, 2.042, 2.124, 0.73
    92.906, 8.57, 154.3, 18.986, 0.468, 16.614, 18.621, 2.166, 2.257, 0.747
    95.95, 10.22, 159.3, 19.999, 0.505, 17.478, 19.607, 2.293, 2.395, 0.765
    99, 11.5, 164.4, 21.045, 0.544, 18.41, 20.585, 2.424, 2.538, 0.78
    101.07, 12.41, 167.8, 22.117, 0.585, 19.278, 21.655, 2.558, 2.683, 0.794
    102.91, 12.44, 170.9, 23.22, 0.627, 20.214, 22.721, 2.696, 2.834, 0.808
    106.4, 12.16, 176.7, 24.35, 0.67, 21.175, 23.816, 2.838, 2.99, 0.82
    107.88, 10.5, 179.1, 25.514, 0.717, 22.162, 24.942, 2.984, 3.151, 0.831
    112.41, 8.65, 186.6, 26.711, 0.77, 23.172, 26.093, 3.133, 3.316, 0.843
    114.82, 7.28, 190.7, 27.94, 0.825, 24.207, 27.274, 3.287, 3.487, 0.853
    118.69, 5.31, 197.1, 29.2, 0.884, 25.27, 28.483, 3.444, 3.662, 0.862
    121.76, 6.691, 202.2, 30.491, 0.944, 26.357, 29.723, 3.605, 3.843, 0.87
    127.6, 6.24, 211.9, 31.813, 1.006, 27.471, 30.993, 3.769, 4.029, 0.877
    126.91, 4.94, 210.7, 33.169, 1.072, 28.61, 32.292, 3.937, 4.22, 0.884
    131.3, 0.0059, 218, 34.582, 1.143, 29.802, 33.644, 4.111, 4.422, 0.891
    132.91, 1.873, 220.7, 35.985, 1.218, 30.97, 34.984, 4.286, 4.62, 0.897
    137.36, 3.5, 228.1, 37.441, 1.293, 32.191, 36.376, 4.467, 4.828, 0.902
    138.92, 6.15, 230.7, 38.925, 1.363, 33.44, 37.799, 4.651, 5.043, 0.907
    140.13, 6.67, 232.7, 40.444, 1.434, 34.717, 39.255, 4.84, 5.262, 0.912
    140.92, 6.769, 234, 41.991, 1.508, 36.023, 40.746, 5.034, 5.489, 0.917
    144.27, 6.96, 239.6, 43.569, 1.575, 37.359, 42.269, 5.431, 5.956, 0.921
    147, 6.782, 244.1, 45.184, 1.651, 38.649, 43.945, 5.636, 6.206, 0.925
    150.35, 7.536, 249.6, 46.835, 1.729, 40.124, 45.4, 5.846, 6.456, 0.929
    152, 5.259, 252.4, 48.52, 1.8, 41.529, 47.027, 6.059, 6.714, 0.932
    157.26, 7.95, 261.1, 50.24, 1.882, 42.983, 48.718, 6.275, 6.979, 0.935
    158.93, 8.272, 263.9, 51.996, 1.967, 44.47, 50.391, 6.495, 7.249, 0.938
    162.51, 8.536, 269.8, 53.789, 2.046, 45.985, 52.178, 6.72, 7.528, 0.941
    164.94, 8.803, 273.9, 55.618, 2.127, 47.528, 53.934, 6.948, 7.81, 0.944
    167.27, 9.051, 277.7, 57.486, 2.212, 49.099, 55.69, 7.181, 8.103, 0.947
    168.94, 9.332, 280.5, 59.39, 2.307, 50.73, 57.576, 7.181, 8.103, 0.949
    173.04, 6.977, 287.3, 61.332, 2.398, 52.36, 59.352, 7.414, 8.401, 0.951
    174.99, 9.842, 290.6, 63.314, 2.492, 54.063, 61.282, 7.654, 8.708, 0.953
    178.5, 13.3, 296.4, 65.351, 2.602, 55.757, 63.209, 7.898, 9.021, 0.955
    180.95, 16.6, 300.5, 67.414, 2.703, 57.524, 65.21, 8.145, 9.341, 0.957
    183.92, 19.3, 305.4, 69.524, 2.818, 59.31, 67.233, 8.396, 9.67, 0.958
    186.2, 21.02, 309.2, 71.676, 2.931, 61.131, 69.298, 8.651, 10.008, 0.959
    190.2, 22.5, 315.8, 73.872, 3.05, 62.991, 71.404, 8.91, 10.354, 0.961
    192.2, 22.42, 319.1, 76.112, 3.172, 64.886, 73.549, 9.173, 10.706, 0.962
    195.09, 21.37, 323.9, 78.395, 3.297, 66.82, 75.736, 9.441, 11.069, 0.963
    197.2, 19.37, 327.4, 80.723, 3.425, 68.779, 77.968, 9.711, 11.439, 0.964
    200.61, 13.546, 333.1, 83.103, 3.562, 70.821, 80.258, 9.987, 11.823, 0.965
    204.39, 11.86, 339.4, 85.528, 3.7, 72.86, 82.558, 10.266, 12.21, 0.966
    207.21, 11.34, 344.1, 88.006, 3.85, 74.957, 84.922, 10.549, 12.611, 0.967
    209, 9.8, 347, 90.527, 3.999, 77.097, 87.335, 10.836, 13.021, 0.968
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.968
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.969
    222, 0.00973, 368.6, 98.417, 4.478, 83.8, 94.877, 11.724, 14.316, 0.969
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.97
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.97
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.971
    232, 11.7, 385.2, 109.649, 5.182, 93.334, 105.592, 12.966, 16.2, 0.971
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.972
    238.07, 19.05, 395.3, 115.603, 5.549, 98.428, 111.289, 13.613, 17.218, 0.972
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0.973
    239.1, 19.7, 397, 121.76, 5.914, 103.653, 117.146, 14.279, 18.278, 0.973 ];

% three-element vectors: L_edge(1:3)  L_jump(1:3)  fluo_L(1:3)
V = [ ...
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0, 0, 0, 1.16, 1.41, 0, 0, 0, 0
    0.063, 0.0496, 0.0492, 1.16, 1.41, 0, 2.9e-05, 0.0012, 0.0012
    0.087, 0.0729, 0.0725, 1.16, 1.41, 0, 2.6e-05, 0.00075, 0.00075
    0.118, 0.0998, 0.0992, 1.16, 1.41, 0, 3e-05, 0.00037, 0.00038
    0.153, 0.136, 0.135, 1.16, 1.41, 0, 3.9e-05, 0.00031, 0.00031
    0.193, 0.1636, 0.1625, 1.16, 1.41, 0, 7.4e-05, 0.00026, 0.00026
    0.238, 0.202, 0.2, 1.16, 1.41, 0, 0.00012, 0.00024, 0.00024
    0.287, 0.2506, 0.2484, 1.16, 1.41, 0, 0.00018, 0.00022, 0.00022
    0.341, 0.2973, 0.2946, 1.16, 1.41, 0, 0.00024, 0.00027, 0.00027
    0.4, 0.3497, 0.3462, 1.16, 1.41, 0, 0.00031, 0.00033, 0.00033
    0.463, 0.4036, 0.3987, 1.16, 1.41, 0, 0.00039, 0.00084, 0.00084
    0.531, 0.4612, 0.4538, 1.16, 1.41, 0, 0.00047, 0.0015, 0.0015
    0.604, 0.5198, 0.5121, 1.16, 1.41, 0, 0.00058, 0.0026, 0.0026
    0.682, 0.5838, 0.5741, 1.16, 1.41, 0, 0.00071, 0.0037, 0.0037
    0.754, 0.6499, 0.6387, 1.16, 1.41, 0, 0.00084, 0.005, 0.005
    0.842, 0.7199, 0.7068, 1.16, 1.41, 0, 0.001, 0.0063, 0.0063
    0.929, 0.7933, 0.7781, 1.16, 1.41, 0, 0.0012, 0.0077, 0.0077
    1.012, 0.872, 0.855, 1.16, 1.41, 2.772, 0.0014, 0.0086, 0.0093
    1.1, 0.952, 0.932, 1.16, 1.41, 2.874, 0.0016, 0.01, 0.011
    1.196, 1.044, 1.021, 1.16, 1.41, 5.684, 0.0018, 0.011, 0.012
    1.302, 1.142, 1.115, 1.16, 1.41, 5.671, 0.0021, 0.012, 0.013
    1.414, 1.249, 1.218, 1.16, 1.41, 5.704, 0.0024, 0.013, 0.015
    1.53, 1.36, 1.325, 1.16, 1.41, 4.875, 0.0028, 0.014, 0.016
    1.653, 1.477, 1.436, 1.16, 1.41, 4.587, 0.0032, 0.016, 0.018
    1.782, 1.596, 1.55, 1.16, 1.41, 4.557, 0.0036, 0.018, 0.02
    1.92, 1.726, 1.675, 1.16, 1.41, 4.17, 0.0041, 0.02, 0.022
    2.065, 1.863, 1.805, 1.16, 1.41, 4.223, 0.0046, 0.022, 0.024
    2.216, 2.007, 1.94, 1.16, 1.41, 3.906, 0.0051, 0.024, 0.026
    2.373, 2.156, 2.08, 1.16, 1.41, 4.036, 0.0059, 0.026, 0.028
    2.532, 2.307, 2.223, 1.16, 1.41, 3.976, 0.0068, 0.028, 0.031
    2.698, 2.465, 2.371, 1.16, 1.41, 3.774, 0.0094, 0.031, 0.034
    2.866, 2.625, 2.52, 1.16, 1.41, 3.675, 0.01, 0.034, 0.037
    3.043, 2.793, 2.677, 1.16, 1.41, 3.591, 0.011, 0.037, 0.04
    3.224, 2.967, 2.838, 1.16, 1.41, 3.431, 0.012, 0.04, 0.043
    3.412, 3.146, 3.003, 1.16, 1.41, 3.721, 0.013, 0.043, 0.046
    3.605, 3.33, 3.173, 1.16, 1.41, 3.402, 0.014, 0.047, 0.049
    3.806, 3.524, 3.351, 1.16, 1.41, 3.223, 0.016, 0.051, 0.052
    4.018, 3.727, 3.537, 1.16, 1.41, 3.249, 0.018, 0.056, 0.056
    4.238, 3.938, 3.73, 1.16, 1.41, 3.255, 0.02, 0.061, 0.06
    4.465, 4.156, 3.929, 1.16, 1.41, 3.06, 0.037, 0.065, 0.064
    4.698, 4.381, 4.132, 1.16, 1.41, 2.939, 0.039, 0.069, 0.069
    4.939, 4.612, 4.341, 1.16, 1.41, 2.979, 0.041, 0.074, 0.074
    5.188, 4.852, 4.557, 1.16, 1.41, 2.856, 0.044, 0.079, 0.079
    5.452, 5.1, 4.781, 1.16, 1.41, 2.879, 0.046, 0.083, 0.085
    5.713, 5.359, 5.012, 1.16, 1.41, 2.847, 0.049, 0.09, 0.091
    5.987, 5.624, 5.247, 1.16, 1.41, 2.839, 0.052, 0.096, 0.097
    6.267, 5.891, 5.483, 1.16, 1.41, 2.717, 0.055, 0.103, 0.104
    6.549, 6.165, 5.724, 1.16, 1.41, 2.737, 0.058, 0.11, 0.111
    6.835, 6.441, 5.965, 1.16, 1.41, 2.695, 0.061, 0.117, 0.118
    7.126, 6.722, 6.208, 1.16, 1.41, 2.662, 0.064, 0.124, 0.125
    7.428, 7.013, 6.46, 1.16, 1.41, 2.702, 0.066, 0.132, 0.132
    7.737, 7.312, 6.717, 1.16, 1.41, 2.68, 0.071, 0.14, 0.139
    8.052, 7.618, 6.977, 1.16, 1.41, 2.723, 0.075, 0.149, 0.147
    8.376, 7.931, 7.243, 1.16, 1.41, 2.701, 0.079, 0.158, 0.155
    8.708, 8.252, 7.515, 1.16, 1.41, 2.713, 0.083, 0.167, 0.164
    9.047, 8.581, 7.79, 1.16, 1.41, 9.047, 0.089, 0.178, 0.174
    9.395, 8.919, 8.071, 1.16, 1.41, 2.863, 0.094, 0.189, 0.182
    9.752, 9.265, 8.358, 1.16, 1.41, 2.933, 0.1, 0.2, 0.192
    10.116, 9.618, 8.648, 1.16, 1.41, 2.758, 0.106, 0.211, 0.201
    10.488, 9.978, 8.943, 1.16, 1.41, 2.573, 0.112, 0.222, 0.21
    10.87, 10.349, 9.244, 1.16, 1.41, 2.62, 0.12, 0.234, 0.22
    11.272, 10.739, 9.561, 1.16, 1.41, 2.415, 0.128, 0.246, 0.231
    11.68, 11.136, 9.881, 1.16, 1.41, 2.6, 0.137, 0.258, 0.243
    12.098, 11.542, 10.204, 1.16, 1.41, 2.617, 0.147, 0.27, 0.255
    12.525, 11.957, 10.534, 1.16, 1.41, 2.675, 0.144, 0.283, 0.268
    12.964, 12.384, 10.871, 1.16, 1.41, 2.529, 0.13, 0.295, 0.281
    13.424, 12.824, 11.215, 1.16, 1.41, 2.387, 0.12, 0.308, 0.294
    13.892, 13.273, 11.564, 1.16, 1.41, 2.632, 0.114, 0.321, 0.306
    14.353, 13.733, 11.918, 1.16, 1.41, 2.439, 0.107, 0.334, 0.32
    14.846, 14.209, 12.284, 1.16, 1.41, 2.4, 0.107, 0.347, 0.333
    15.344, 14.698, 12.657, 1.16, 1.41, 2.498, 0.107, 0.36, 0.347
    15.86, 15.198, 13.035, 1.16, 1.41, 2.466, 0.112, 0.373, 0.36
    16.385, 15.708, 13.418, 1.16, 1.41, 2.338, 0.117, 0.387, 0.373
    0, 0, 0, 1.16, 1.41, 0, 0.122, 0.401, 0.386
    0, 0, 0, 1.16, 1.41, 0, 0.128, 0.415, 0.399
    18.055, 17.334, 14.612, 1.16, 1.41, 2.344, 0.134, 0.429, 0.411
    0, 0, 0, 1.16, 1.41, 0, 0.139, 0.443, 0.424
    0, 0, 0, 1.16, 1.41, 0, 0.146, 0.456, 0.437
    0, 0, 0, 1.16, 1.41, 0, 0.153, 0.468, 0.45
    20.47, 19.692, 16.3, 1.16, 1.41, 2.388, 0.161, 0.479, 0.463
    0, 0, 0, 1.16, 1.41, 0, 0.162, 0.472, 0.476
    21.756, 20.947, 17.167, 1.16, 1.41, 2.292, 0.176, 0.467, 0.489
    0, 0, 0, 1.16, 1.41, 0, 0.187, 0.466, 0.502
    23.095, 22.263, 18.053, 1.16, 1.41, 2.251, 0.205, 0.464, 0.514 ];

% four-element vectors, in this order:
%   xsec_coh  xsec_nch  post_K  post_L  post_M  post_N
P = [ ...
    -0.119075, -0.937086, -0.200538, 0.0106587, -2.15772, 1.32685, -0.30562, 0.0185025, 2.44964, -3.34953, -0.047137, 0.00709962, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    1.04768, -0.0851805, -0.403527, 0.0269398, -2.56357, 2.02536, -0.44871, 0.0279691, 6.06488, -3.29055, -0.107256, 0.0144465, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    1.34366, 0.181557, -0.423981, 0.026619, -1.0874, 1.03368, -0.190377, 0.00779955, 7.7537, -2.81801, -0.241378, 0.0262542, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    2.0086, -0.046192, -0.337018, 0.0186939, -0.690079, 0.946448, -0.171142, 0.00651413, 9.04511, -2.83487, -0.210021, 0.0229526, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    2.62862, -0.207916, -0.286283, 0.0144966, -0.791177, 1.21611, -0.239087, 0.0117686, 9.95057, -2.74174, -0.215138, 0.0227845, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    3.10861, -0.26058, -0.271974, 0.0135181, -0.982878, 1.46693, -0.293743, 0.0156005, 10.6879, -2.714, -0.20053, 0.0207248, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    3.4776, -0.215762, -0.288874, 0.0115131, -1.23693, 1.7451, -0.35466, 0.0198705, 11.2765, -2.654, -0.200445, 0.0200765, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    3.77239, -0.148539, -0.307124, 0.0167303, -1.73679, 2.17686, -0.44905, 0.0264733, 11.713, -2.57229, -0.205893, 0.0199244, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    4.00716, -0.0560908, -0.332017, 0.0187934, -1.8757, 2.32016, -0.475412, 0.028068, 12.0963, -2.44148, -0.234461, 0.0219954, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    4.20151, 0.0416247, -0.356754, 0.0207585, -1.7551, 2.24226, -0.44764, 0.0255801, 12.4485, -2.45819, -0.212591, 0.0196489, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    4.26374, 0.134662, -0.37008, 0.0214467, -0.967717, 1.61794, -0.287191, 0.0131526, 12.6777, -2.24521, -0.274873, 0.025027, 10.2355, -2.55905, -0.119524, 0, 0, 0, 0, 0, 0, 0, 0, 0
    4.39404, 0.137858, -0.35954, 0.020238, -0.571611, 1.35498, -0.222491, 0.00830141, 12.8793, -2.12574, -0.299392, 0.0267643, 10.5973, -2.89818, 0.234506, 0, 0, 0, 0, 0, 0, 0, 0, 0
    4.15995, 0.140549, -0.352441, 0.0193692, -0.439322, 1.30867, -0.211648, 0.0075421, 13.1738, -2.18203, -0.25896, 0.022284, 10.8711, -2.7786, 0.175853, 0, 0, 0, 0, 0, 0, 0, 0, 0
    4.64678, 0.16278, -0.358563, 0.0196926, -0.414971, 1.34868, -0.222315, 0.00841959, 13.2682, -1.98174, -0.31695, 0.0273928, 11.2237, -2.73694, 0.127557, 0, 0, 0, 0, 0, 0, 0, 0, 0
    4.78525, 0.168708, -0.360383, 0.0197155, -0.476903, 1.46032, -0.251331, 0.0107202, 13.3735, -1.86342, -0.33944, 0.0288858, 11.5508, -2.922, 0.254262, 0, 0, 0, 0, 0, 0, 0, 0, 0
    4.92707, 0.165746, -0.359424, 0.0195505, -0.656419, 1.65408, -0.298623, 0.0142979, 13.7394, -2.04786, -0.273259, 0.0229976, 11.8181, -2.64618, -0.0968049, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.07222, 0.149127, -0.352858, 0.0189439, -0.718627, 1.74294, -0.319429, 0.0158429, 13.6188, -1.71937, -0.354154, 0.0290841, 12.0031, -2.41694, -0.240897, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.21079, 0.135618, -0.347214, 0.0184333, -0.682105, 1.74279, -0.317646, 0.0156467, 13.9491, -1.82276, -0.328827, 0.0274382, 12.296, -2.63279, -0.07366, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.25587, 0.18804, -0.359623, 0.0193085, -0.344007, 1.49236, -0.254135, 0.0107684, 13.7976, -1.54015, -0.394528, 0.0323561, 12.4878, -2.53656, -0.104892, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.32375, 0.206685, -0.361664, 0.0193328, -0.098242, 1.32829, -0.213747, 0.00773065, 14.295, -1.88644, -0.283647, 0.0226263, 12.7044, -2.55011, -0.0943195, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.43942, 0.200174, -0.359064, 0.0191027, -0.159831, 1.39055, -0.225849, 0.00851954, 13.9664, -1.40872, -0.414365, 0.0334355, 12.8949, -2.40609, -0.177791, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.55039, 0.197697, -0.357694, 0.0189866, -0.230573, 1.45848, -0.23916, 0.00938528, 14.3506, -1.66322, -0.331539, 0.0262065, 13.1075, -2.53576, -0.0957177, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.65514, 0.199533, -0.357487, 0.0189691, -0.308103, 1.52879, -0.252768, 0.0102571, 14.7601, -1.88867, -0.271861, 0.0215792, 13.2514, -2.49765, -0.106383, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.77399, 0.203858, -0.359699, 0.0192225, -0.387641, 1.59727, -0.26624, 0.0111523, 14.8019, -1.8243, -0.279116, 0.0217324, 13.4236, -2.51532, -0.101999, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.84604, 0.213814, -0.359718, 0.0191459, -0.247059, 1.49722, -0.238781, 0.00893208, 14.8965, -1.79872, -0.283664, 0.0222095, 13.5761, -2.49761, -0.105493, 0, 0, 0, 0, 0, 0, 0, 0, 0
    5.93292, 0.225038, -0.361748, 0.0193024, -0.342379, 1.57245, -0.253198, 0.00985822, 14.3456, -1.23491, -0.423491, 0.0321661, 13.6696, -2.39195, -0.13768, 0, 0, 0, 0, 0, 0, 0, 0, 0
    6.01478, 0.237959, -0.364056, 0.0194754, -0.428804, 1.64129, -0.266013, 0.0106512, 14.7047, -1.38933, -0.386631, 0.0303286, 13.8699, -2.50669, -0.0869945, 0, 0, 0, 0, 0, 0, 0, 0, 0
    6.09024, 0.22277, -0.366568, 0.0196586, -0.50436, 1.7004, -0.276443, 0.0112628, 14.2388, -0.967736, -0.47807, 0.0366138, 13.9848, -2.4808, -0.0888115, 0, 0, 0, 0, 0, 0, 0, 0, 0
    6.17739, 0.273123, -0.37236, 0.0201638, -0.57021, 1.75042, -0.284555, 0.011693, 14.5808, -1.18375, -0.41385, 0.0312088, 14.2439, -2.58677, -0.0667398, 0, 0, 0, 0, 0, 0, 0, 0, 0
    6.23402, 0.284312, -0.372143, 0.0200525, -0.420535, 1.634, -0.253646, 0.00927233, 14.4118, -0.933083, -0.477357, 0.0362829, 14.3221, -2.62384, -0.0264926, 0, 12.0597, -1.10258, 0, 0, 0, 0, 0, 0
    6.28298, 0.291334, -0.369391, 0.0197029, -0.358218, 1.6005, -0.244908, 0.00861898, 13.6182, -0.318459, -0.611348, 0.0458138, 14.4792, -2.54469, -0.0757204, 0, 12.2646, -2.68965, 0, 0, 0, 0, 0, 0
    6.33896, 0.291512, -0.365643, 0.0192896, -0.334383, 1.60327, -0.245555, 0.00871239, 13.9288, -0.479613, -0.572897, 0.0431277, 14.6813, -2.69285, -0.0208355, 0, 12.4133, -2.53085, 0, 0, 0, 0, 0, 0
    6.3975, 0.288866, -0.361747, 0.0188788, -0.339189, 1.62535, -0.250783, 0.00909103, 13.4722, -0.0773513, -0.660456, 0.0492177, 14.6431, -2.48397, -0.079618, 0, 12.5392, -2.4138, 0, 0, 0, 0, 0, 0
    6.45637, 0.286737, -0.358794, 0.0185618, -0.432927, 1.72833, -0.277138, 0.0111735, 13.0756, 0.183235, -0.694264, 0.050228, 14.7048, -2.38853, -0.105877, 0, 12.6773, -2.3975, 0, 0, 0, 0, 0, 0
    6.51444, 0.286324, -0.357027, 0.0183557, -0.448001, 1.76082, -0.285099, 0.0117865, 13.2273, 0.13713, -0.683203, 0.0495424, 14.8136, -2.42347, -0.091459, 0, 12.7612, -2.3773, 0, 0, 0, 0, 0, 0
    6.57113, 0.287711, -0.356311, 0.018247, -0.39181, 1.7301, -0.276824, 0.011128, 13.5927, -0.0305214, -0.65134, 0.0477616, 14.919, -2.42418, -0.0876447, 0, 12.8898, -2.26021, 0, 0, 0, 0, 0, 0
    6.5975, 0.302389, -0.356755, 0.0181706, -0.128039, 1.53044, -0.227403, 0.00739033, 13.0204, 0.382736, -0.732427, 0.0529874, 14.9985, -2.39108, -0.0959473, 0, 13.0286, -2.38693, 0, 0, 0, 0, 0, 0
    6.62203, 0.324559, -0.361651, 0.01848, 0.0799161, 1.38397, -0.192225, 0.00478611, 13.5888, 0.00220194, -0.63894, 0.046007, 15.0114, -2.28169, -0.126485, 0, 13.1565, -2.36655, 0, 0, 0, 0, 0, 0
    6.67096, 0.32507, -0.360613, 0.0183326, 0.0629057, 1.41577, -0.199713, 0.00533312, 13.4674, 0.191023, -0.686616, 0.0497356, 15.1822, -2.38946, -0.0881174, 0, 13.2775, -2.43174, 0, 0, 0, 0, 0, 0
    6.72275, 0.323964, -0.359463, 0.018189, 0.0366697, 1.45207, -0.208122, 0.00595139, 12.7538, 0.697409, -0.789307, 0.0564531, 15.2906, -2.38703, -0.0912292, 0, 13.4508, -2.50201, 0, 0, 0, 0, 0, 0
    6.79013, 0.311282, -0.355233, 0.0178231, 0.000202289, 1.49347, -0.217419, 0.00662245, 13.3843, 0.281028, -0.686607, 0.0486607, 15.2088, -2.20278, -0.136759, 0, 13.5434, -2.50135, 0, 0, 0, 0, 0, 0
    6.846, 0.30279, -0.351131, 0.0174403, -0.056286, 1.55778, -0.233341, 0.00785506, 13.9853, -0.117426, -0.591094, 0.0417843, 15.3494, -2.2664, -0.116881, 0, 13.6568, -2.4848, 0, 0, 0, 0, 0, 0
    6.87599, 0.326165, -0.358969, 0.0180482, 0.0757616, 1.4495, -0.20489, 0.00564745, 12.8214, 0.751993, -0.787006, 0.0558668, 15.5086, -2.33733, -0.0987857, 0, 13.7498, -2.4473, 0, 0, 0, 0, 0, 0
    6.93136, 0.334794, -0.363497, 0.0184429, -0.0424981, 1.54639, -0.22647, 0.00718375, 12.6658, 0.88502, -0.811144, 0.0573759, 15.4734, -2.2308, -0.119454, 0, 13.8782, -2.48066, 0, 0, 0, 0, 0, 0
    6.97547, 0.346394, -0.367794, 0.0187885, -0.160399, 1.64861, -0.250238, 0.00893818, 12.176, 1.19682, -0.866697, 0.0606931, 15.5757, -2.24976, -0.113377, 0, 14.0312, -2.61303, 0, 0, 0, 0, 0, 0
    7.03216, 0.349838, -0.370099, 0.0189983, -0.267564, 1.7374, -0.269883, 0.0103248, 13.9389, 0.164528, -0.662117, 0.0476289, 15.5649, -2.17229, -0.127652, 0, 14.1392, -2.57206, 0, 0, 0, 0, 0, 0
    7.06446, 0.363456, -0.373597, 0.0192478, -0.166475, 1.65794, -0.24874, 0.00866218, 13.3926, 0.44138, -0.693711, 0.0482085, 15.6869, -2.22636, -0.112223, 0, 14.1673, -2.48078, 0, 0, 0, 0, 0, 0
    7.09856, 0.372199, -0.375345, 0.0193481, -0.0516701, 1.57426, -0.227646, 0.0070565, 12.5254, 1.07714, -0.831424, 0.057912, 15.9668, -2.38363, -0.0801104, 0, 14.3497, -2.52756, 0, 0, 0, 0, 0, 0
    7.12708, 0.382082, -0.376855, 0.0194151, -0.00817283, 1.55865, -0.224492, 0.00685776, 11.8198, 1.45768, -0.888529, 0.0605982, 16.2101, -2.51838, -0.0540061, 0, 14.4115, -2.49401, 0, 0, 0, 0, 0, 0
    7.16085, 0.385512, -0.376481, 0.0193305, 0.0142151, 1.55754, -0.224736, 0.00691395, 13.0323, 0.790788, -0.762349, 0.0527872, 15.8638, -2.19019, -0.113539, 0, 14.5572, -2.56792, 0, 0, 0, 0, 0, 0
    7.19665, 0.385543, -0.375054, 0.0191608, 0.0156362, 1.57175, -0.228753, 0.00726386, 9.0699, 3.28791, -1.26203, 0.085347, 15.7557, -2.0446, -0.140745, 0, 14.6268, -2.55562, 0, 0, 0, 0, 0, 0
    7.2346, 0.382493, -0.372715, 0.0189194, -0.0407579, 1.64267, -0.24789, 0.00880567, 11.6656, 1.71052, -0.948281, 0.0653213, 16.1087, -2.27876, -0.0929405, 0, 14.7125, -2.54324, 0, 0, 0, 0, 0, 0
    7.27415, 0.377223, -0.369728, 0.018628, -0.040442, 1.65596, -0.251067, 0.00904874, 12.1075, 1.43635, -0.882038, 0.0603575, 16.4086, -2.48214, -0.0507179, 0, 14.7496, -2.48179, 0, 0, 0, 0, 0, 0
    7.31469, 0.370315, -0.36628, 0.0183025, -0.00282407, 1.64039, -0.247642, 0.00882144, 11.0857, 2.08357, -1.01209, 0.069031, 16.3098, -2.31679, -0.0854498, 0, 14.7603, -2.45068, 0, 0, 0, 0, 0, 0
    7.3349, 0.376825, -0.365713, 0.0181843, 0.184861, 1.5003, -0.213333, 0.00624264, 11.375, 1.94161, -0.983232, 0.0671986, 16.5418, -2.46363, -0.0542849, 0, 14.9713, -2.53145, 0, 0, 0, 0, 0, 0
    7.35812, 0.379361, -0.364099, 0.0179817, 0.344376, 1.38742, -0.186356, 0.00424917, 10.225, 2.67835, -1.12648, 0.0762669, 16.6217, -2.48972, -0.0449623, 0, 15.0844, -2.56341, 0, 0, 0, 0, 0, 0
    7.39532, 0.369895, -0.359376, 0.0175406, 0.409104, 1.33075, -0.170883, 0.00304111, 10.978, 2.23814, -1.03549, 0.0702339, 16.3134, -2.20156, -0.0980569, 0, 15.1863, -2.58287, 0, 0, 0, 0, 0, 0
    7.44255, 0.371328, -0.359642, 0.0175852, 0.439881, 1.30925, -0.164548, 0.00252641, 10.2725, 2.74562, -1.14174, 0.0774162, 16.5862, -2.36288, -0.0654708, 0, 15.2693, -2.58174, 0, 0, 0, 0, 0, 0
    7.48347, 0.368431, -0.357689, 0.0174099, 0.449124, 1.30351, -0.161841, 0.00227394, 11.0156, 2.22056, -1.02216, 0.0690465, 16.7118, -2.40326, -0.0612619, 0, 15.3379, -2.57086, 0, 0, 0, 0, 0, 0
    7.52334, 0.366462, -0.356048, 0.017262, 0.437283, 1.3137, -0.162866, 0.00229377, 11.7632, 1.79481, -0.936661, 0.0635332, 16.5964, -2.26007, -0.0872426, 0, 15.4335, -2.59006, 0, 0, 0, 0, 0, 0
    7.56222, 0.365055, -0.354511, 0.0171214, 0.405823, 1.33837, -0.167229, 0.0025557, 11.3864, 2.05593, -0.98818, 0.0669106, 16.8337, -2.3881, -0.0645041, 0, 15.5131, -2.59623, 0, 0, 15.5131, -2.59623, 0, 0
    7.6002, 0.364134, -0.353086, 0.0169894, 0.355383, 1.37733, -0.174941, 0.00306213, 11.9223, 1.79546, -0.942902, 0.0644202, 16.8725, -2.39051, -0.060108, 0, 15.6, -2.61328, 0, 0, 15.6006, -2.61328, 0, 0
    7.63711, 0.363957, -0.3519, 0.0168783, 0.280316, 1.44016, -0.188641, 0.00401226, 11.6168, 1.97533, -0.970901, 0.0658459, 17.0692, -2.48046, -0.0447055, 0, 15.7063, -2.63481, 0, 0, 15.7063, -2.63481, 0, 0
    7.66938, 0.359752, -0.348899, 0.016589, 0.273133, 1.43842, -0.186137, 0.0037524, 9.91968, 3.03111, -1.1752, 0.078675, 17.1159, -2.47838, -0.0437107, 0, 15.7159, -2.60843, 0, 0, 15.7159, -2.60843, 0, 0
    7.70798, 0.365345, -0.350031, 0.0166927, 0.257539, 1.45064, -0.187591, 0.00379932, 11.3818, 2.14447, -0.999222, 0.0675569, 17.1499, -2.45507, -0.047137, 0, 15.8412, -2.6404, 0, 0, 15.8415, -2.6404, 0, 0
    7.74188, 0.367107, -0.349433, 0.0166273, 0.242685, 1.46266, -0.189102, 0.00385628, 11.4845, 2.10451, -0.98987, 0.0669382, 17.3446, -2.54821, -0.0317606, 0, 15.9225, -2.65289, 0, 0, 15.9225, -2.65289, 0, 0
    7.7747, 0.369722, -0.349132, 0.0165862, 0.228493, 1.47438, -0.190559, 0.00390903, 8.75203, 3.71822, -1.29273, 0.0855026, 17.6583, -2.72523, -0.000819409, 0, 16.014, -2.67903, 0, 0, 16.014, -2.67903, 0, 0
    7.80643, 0.373226, -0.349147, 0.016571, 0.225233, 1.48545, -0.191908, 0.00395645, 12.0195, 1.84815, -0.939582, 0.0638106, 17.7988, -2.74671, -0.0028758, 0, 16.0672, -2.6758, 0, 0, 16.0672, -2.6758, 0, 0
    7.83711, 0.377547, -0.349441, 0.016578, 0.202656, 1.49625, -0.193234, 0.00400233, 12.5613, 1.57523, -0.890467, 0.0609779, 17.425, -2.51103, -0.032945, 0, 16.1269, -2.67886, 0, 0, 16.1269, -2.67886, 0, 0
    7.86662, 0.382933, -0.350126, 0.0166173, 0.202248, 1.48804, -0.189143, 0.00362264, 7.42791, 4.28955, -1.35167, 0.0866136, 16.9795, -2.22577, -0.0732557, 0, 16.1794, -2.67715, 0, 0, 13.9111, -2.4038, 0, 0
    7.89137, 0.386034, -0.349756, 0.016548, 0.197176, 1.50264, -0.192474, 0.00385751, 12.6387, 1.55476, -0.881094, 0.0602036, 17.2638, -2.37189, -0.0495994, 0, 16.2289, -2.67128, 0, 0, 13.9813, -2.40841, 0, 0
    7.91803, 0.387021, -0.348881, 0.0164406, 0.199469, 1.50233, -0.191385, 0.00374011, 7.5816, 4.47037, -1.42808, 0.0939044, 16.4329, -1.82851, -0.132268, 0, 16.2758, -2.6622, 0, 0, 14.0548, -2.42829, 0, 0
    7.94534, 0.387299, -0.347926, 0.016323, 0.196871, 1.50623, -0.191396, 0.00370889, 8.65271, 3.73117, -1.26359, 0.0823539, 17.2411, -2.30313, -0.0591006, 0, 16.3038, -2.66148, 0, 0, 14.1313, -2.47214, 0, 0
    7.97266, 0.387704, -0.347155, 0.0162372, 0.191015, 1.5124, -0.19122, 0.0037145, 7.57541, 4.28874, -1.34998, 0.08652, 17.2533, -2.23874, -0.0727338, 0, 16.2613, -2.60672, 0, 0, 14.1536, -2.3258, 0, 0
    7.9994, 0.388739, -0.346726, 0.0161751, 0.189644, 1.50867, -0.18957, 0.00349584, 1.36944, 7.79444, -1.99822, 0.126225, 17.875, -2.61051, -0.0136093, 0, 16.3564, -2.62453, 0, 0, 14.2392, -2.35326, 0, 0
    8.02574, 0.390458, -0.346658, 0.0161455, 0.116448, 1.57615, -0.205332, 0.00466731, 13.7534, 1.02122, -0.777126, 0.0538811, 17.3525, -2.2855, -0.0588047, 0, 16.4233, -2.63163, 0, 0, 14.2795, -2.21971, 0, 0
    8.0515, 0.393143, -0.347052, 0.016157, 0.0719908, 1.61204, -0.213186, 0.00520497, 12.5506, 1.6309, -0.875676, 0.0592011, 16.527, -1.76315, -0.135232, 0, 16.5144, -2.64832, 0, 0, 14.3422, -2.40183, 0, 0
    8.08084, 0.39579, -0.348032, 0.0162345, 0.0420186, 1.63611, -0.217964, 0.0055267, 12.7882, 1.63605, -0.898523, 0.061855, 17.3636, -2.2112, -0.0730934, 0, 16.7024, -2.71631, 0, 0, 14.3785, -2.34834, 0, 0
    8.10524, 0.400576, -0.34934, 0.0163264, 0.0156916, 1.65406, -0.220982, 0.00570751, 4.96352, 5.79212, -1.61842, 0.102911, 17.424, -2.23911, -0.066372, 0, 16.4734, -2.57834, 0, 0, 14.4398, -2.32838, 0, 0
    8.12542, 0.405858, -0.350329, 0.0163772, 0.114587, 1.58076, -0.202968, 0.00435692, 19.7594, -1.9799, -0.276981, 0.0268856, 17.1857, -2.0847, -0.0853294, 0, 16.5903, -2.6067, 0, 0, 14.5195, -2.33016, 0, 0
    8.14399, 0.408692, -0.349802, 0.016288, 0.147052, 1.56695, -0.200347, 0.00420901, 15.2879, 0.273664, -0.63889, 0.0457495, 17.7379, -2.37745, -0.0433223, 0, 16.6564, -2.61593, 0, 0, 14.5473, -2.26773, 0, 0
    8.15996, 0.418031, -0.35233, 0.016466, 0.182167, 1.54661, -0.195973, 0.00390772, 8.63374, 3.694, -1.21312, 0.0774601, 17.7963, -2.37691, -0.0455883, 0, 16.7131, -2.61538, 0, 0, 14.5771, -2.25279, 0, 0
    8.17489, 0.427916, -0.355068, 0.0166601, 0.18986, 1.56125, -0.200932, 0.00436768, 9.44293, 3.44965, -1.19886, 0.0783484, 17.5348, -2.23353, -0.0596161, 0, 16.7078, -2.58648, 0, 0, 14.6832, -2.3094, 0, 0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    8.22553, 0.451478, -0.362056, 0.0171556, 0.196619, 1.6008, -0.2138, 0.00551717, 15.1782, 0.34902, -0.637638, 0.0451377, 17.5028, -2.14876, -0.0724638, 0, 16.9, -2.60945, 0, 0, 14.7243, -2.12905, 0, 0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    8.27843, 0.479056, -0.367657, 0.0174621, 0.17089, 1.65561, -0.229702, 0.00692516, 13.4336, 1.34805, -0.81328, 0.0555664, 18.5481, -2.61281, -0.00790574, 0, 17.0483, -2.58569, 0, 0, 14.773, -1.91192, 0, 0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    8.3301, 0.478314, -0.36725, 0.0174129, 0.108277, 1.74158, -0.254104, 0.00895056, 13.7951, 1.23983, -0.801545, 0.0553596, 17.5258, -207237, -0.0723932, 0, 17.0353, -2.56903, 0, 0, 14.9036, -2.12148, 0, 0
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    8.38174, 0.477085, -0.366556, 0.0173422, 0.0388791, 1.82229, -0.276099, 0.0107392, 18.2787, -1.17371, -0.368344, 0.0298738, 17.5519, -2.02162, -0.082294, 0, 17.2953, -2.62164, 0, 0, 14.8535, -1.87733, 0, 0 ];

% element.mat stores post_L(2) = -207237 for uranium.  Thorium and
% plutonium sit near -2.1, so that is almost certainly a misplaced
% decimal point in the 1969 table, and it makes the U L-edge
% photo-absorption nonsense.  The value is reproduced above exactly as
% stored; uncomment the next line to use the value it was probably
% meant to be.
% P(92,14) = -2.07237;

for Z = 94:-1:1                      % backwards: allocates in one go
    el(Z).symbol   = sym{Z};
    el(Z).weight   = A(Z,1);
    el(Z).density  = A(Z,2);
    el(Z).m2a      = A(Z,3);
    el(Z).K_edge   = A(Z,4);
    el(Z).L_edge   = V(Z,1:3);
    el(Z).M_edge   = A(Z,5);
    el(Z).K_alpha1 = A(Z,6);
    el(Z).K_beta1  = A(Z,7);
    el(Z).L_alpha1 = A(Z,8);
    el(Z).L_beta1  = A(Z,9);
    el(Z).L_jump   = V(Z,4:6);
    el(Z).fluo_K   = A(Z,10);
    el(Z).fluo_L   = V(Z,7:9);
    el(Z).xsec_coh = P(Z,1:4);
    el(Z).xsec_nch = P(Z,5:8);
    el(Z).post_K   = P(Z,9:12);
    el(Z).post_L   = P(Z,13:16);
    el(Z).post_M   = P(Z,17:20);
    el(Z).post_N   = P(Z,21:24);
end

EL = el;
end

% ======================================================================
%      EMBEDDED TOOLBOX FUNCTIONS   (were McMaster.m, zymbol.m,
%      unify.m and legpoly.m; local copies, used in preference to
%      anything of the same name on the MATLAB path)
% ======================================================================
function [mu_tot, mu_xtr, mu_ph, mu_coh, mu_nch, info] = McMaster(Z, E, unit)
%MCMASTER  Cross-sections from the McMaster tables, May 1969 edition.
%
%   mu = McMaster(Z,E)     total cross-section, [cm^2/g]
%   mu = McMaster(Z,E,'b') the same in [barns/atom]
%
%   Z:  atomic number or symbol (Z_max = 94)
%   E:  photon energies, in keV, ascending
%
%   Tsu-Chien Weng, 07/01/98.  Copyleft, the Penner-Hahn research group.
%   10/24/00: correction on L-edge jumps.
%   2005-03-11: clean up, element.mat remade with make_element.m
%   Embedded in norm2poly_gui: the table now comes from
%   n2p_elementTable, and the three "no data" cases raise errors instead
%   of printing and returning nothing.

Z_max = 94;

if nargin < 1, error('McMaster:noZ','Missing arguments: atomic number.'); end
if ischar(Z), Z = zymbol(Z); end
if nargin < 2, error('McMaster:noE','2nd argument: energies in keV.'); end
if nargin < 3, unit = 'c'; end

element = n2p_elementTable();

if isempty(Z) || numel(Z) ~= 1 || Z < 1 || Z ~= fix(Z)
    error('McMaster:badZ','The atomic number should be a positive integer.')
elseif Z > Z_max
    error('McMaster:heavy','Your atom is too heavy (Z_max = %d).',Z_max)
elseif any(Z == [84 85 87 88 89 91 93])
    % Po, At, Fr, Ra, Ac, Pa, Np
    error('McMaster:noData','McMaster did not tabulate element Z = %d.',Z)
end

E = E(:);

% 1. photo-absorption, in [barns/atom]
idx_K  = ( E >= element(Z).K_edge );
idx_L1 = ( E >= element(Z).L_edge(1) & E <  element(Z).K_edge );
idx_L2 = ( E >= element(Z).L_edge(2) & E <  element(Z).L_edge(1) );
idx_L3 = ( E >= element(Z).L_edge(3) & E <  element(Z).L_edge(2) );
idx_M  = ( E >= element(Z).M_edge    & E <  element(Z).L_edge(3) );

if sum(idx_K|idx_L1|idx_L2|idx_L3|idx_M) ~= numel(E)
    error('McMaster:range', ...
        'Energies outside the tabulated range (below the M edge, %.4f keV).', ...
        element(Z).M_edge);
end

mu_ph = exp([ ...
    polyval(element(Z).post_M(4:-1:1),log(E(idx_M))); ...
    polyval(element(Z).post_L(4:-1:1),log(E(idx_L3|idx_L2|idx_L1))); ...
    polyval(element(Z).post_K(4:-1:1),log(E(idx_K)))]);
mu_ph(idx_L2) = mu_ph(idx_L2)/element(Z).L_jump(1);
mu_ph(idx_L3) = mu_ph(idx_L3)/element(Z).L_jump(1)/element(Z).L_jump(2);

% 2. coherent scattering
mu_coh = exp(polyval(element(Z).xsec_coh(4:-1:1),log(E)));

% 3. incoherent scattering
mu_nch = exp(polyval(element(Z).xsec_nch(4:-1:1),log(E)));

mu_xtr = mu_coh + mu_nch;
mu_tot = mu_ph + mu_xtr;

switch unit
    case {'b','B'}
        % atomic absorption coefficient, [barns/atom]: nothing to do
    otherwise
        % mass absorption coefficient, [cm^2/g]
        mu_ph  = mu_ph /element(Z).m2a;
        mu_coh = mu_coh/element(Z).m2a;
        mu_nch = mu_nch/element(Z).m2a;
        mu_xtr = mu_xtr/element(Z).m2a;
        mu_tot = mu_tot/element(Z).m2a;
end

info = element(Z);
end

% ----------------------------------------------------------------------
function out = zymbol(in)
%ZYMBOL     Convert between atomic number and atomic symbol
%
%   Examples:   zymbol('Mn') gives 25.
%               zymbol(82)   gives 'Pb'.
%
%   Tsu-Chien Weng, 2002-10-03.  Copyleft, the Penner-Hahn research group.
%   Embedded in norm2poly_gui: findstr replaced by strfind.

row1 = 'H He';
row2 = 'LiBeB C N O F Ne';
row3 = 'NaMgAlSiP S ClAr';
row4 = 'K CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKr';
row5 = 'RbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI Xe';
row6 = 'CsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRn';
row7 = 'FrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMt';

periodicTable = [row1 row2 row3 row4 row5 row6 row7];
lastZ = length(periodicTable)/2;

if ischar(in)
    if length(in) == 1, in = [in ' ']; end
    in = lower(in);
    in(1) = upper(in(1));
    hit = strfind(periodicTable,in);
    hit = hit(mod(hit,2) == 1);         % symbols start on odd characters
    out = (hit-1)/2 + 1;
elseif isnumeric(in)
    if any(in > lastZ), error('zymbol:tooLarge','Check the atomic number.  It is too large.'), end
    out = deblank(periodicTable([2*in(:)-1 2*in(:)]));
else
    error('zymbol:badInput','Unrecognized argument format.')
end
end

% ----------------------------------------------------------------------
function [a, mid, span] = unify(A,scale)
%UNIFY  transform [min(A), max(A)] to [-scale, +scale] columnwise.
%       default range: [-1,1].  A constant column is copied through
%       unchanged, to avoid NaNs.
%
%   Tsu-Chien Weng, 07/10/98.  Copyleft, the Penner-Hahn research group.
%   Embedded in norm2poly_gui: the constant-column test now looks at the
%   column being processed rather than the whole span vector.

[m,n] = size(A);
if m == 1, n = 1; A = A'; end

for i = 1:n
    mid(i)  = (max(A(:,i))+min(A(:,i)))/2; %#ok<AGROW>
    span(i) = (max(A(:,i))-min(A(:,i)))/2; %#ok<AGROW>

    if span(i) == 0
        a(:,i) = A(:,i);                   %#ok<AGROW> don't touch constant
    else
        a(:,i) = (A(:,i)-mid(i))/span(i);  %#ok<AGROW>
    end
end

if nargin == 2
    a = a*scale;
end
end

% ----------------------------------------------------------------------
function P_n = legpoly(n)
%LEGPOLY  Legendre polynomial coefficients.
%   (n+1)*P_n+1(x) - (2n+1)*x*P_n(x) + n*P_n-1(x) = 0
%
%   Tsu-Chien Weng, 07/16/98.  Copyleft, the Penner-Hahn research group.

if n < 0, error('legpoly:negative','Not support negative integer!'); end
if (n-floor(n)) ~= 0, error('legpoly:notInteger','Input argument should be an integer!'); end

if n == 0
    P_n = 1;
elseif n == 1
    P_n = [1 0];
else
    P_n2 = 1;
    P_n1 = [1 0];
    for idx = 2:n
        P_n  = 1/idx*((2*idx-1)*[P_n1 0]-(idx-1)*[0 0 P_n2]);
        P_n2 = P_n1;
        P_n1 = P_n;
    end
end
end
