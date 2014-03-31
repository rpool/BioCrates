% optional params:
% nodelabels  => cell array of strings of node labels, def: A,B,...
% nodecolors  => Nx3 matrix of color values in RGB MATLAB format
% nodegroups  => Nx1 vector of node groups, equal values will be grouped in yED, group=0 means no grouping
% coordinates => document!
% directed    => T/F whether graph is directed, def: F
% shrink      => remove unconnected nodes, def: T
% useweights  => represent weights as line widths T/F, def: F
% scaleweight => whether edges are automatically scaled to [0 7] or not, def: T
% customscale => [from to] custom scale of matrix values to min/max edge widths
% edgelabels  => sparse_cell array of edge labels
% edgecolors  => Kx3 sparse matrix in RGB MATLAB format, missing ones are black
% edgecolorA  => NxN sparse matrix attributing each edge to a color
% nodeshapes  => Nx1 vector of node shapes
%   possible node shapes are: 
%   0=ellipse, 1=diamond, 2=rectangle, 3=triangle
%   
%
% by JK, JB
% Version: 2012-06-07
function writeYED(outfile,varargin)

% if the first argument is not a structure, it's A
if ~isstruct(varargin{1})
    A=varargin{1};
    varargin=varargin(2:end);
else
    A=varargin{1}.A;
end

% pre-defined shapes
allshapes={'ellipse','diamond','rectangle','triangle'};

% parse parameters
p = inputParser;
p.StructExpand=true; % if we allow parameter  structure expanding
p.addRequired('outfile', @ischar);
p.addRequired('A', @isnumeric);
p.addParamValue('nodelabels', '', @iscellstr);
p.addParamValue('nodecolors',[],@isnumeric);
p.addParamValue('nodegroups',[],@isnumeric);
p.addParamValue('directed', 0, @islogical);
p.addParamValue('shrink', true, @islogical);
p.addParamValue('useweights', false, @islogical);
p.addParamValue('scaleweights', true, @islogical);
p.addParamValue('customscale', [], @isnumeric);
p.addParamValue('edgecolors',[],@isnumeric);
p.addParamValue('edgecolorA',[],@isnumeric);
p.addParamValue('coordinates',[],@isnumeric);
p.addParamValue('edgelabels',[],@(x)isa(x,'sparse_cell'));
p.addParamValue('nodeshapes',[],@isnumeric);
p.addParamValue('nodeshapelist',[],@iscellstr);
p.parse(outfile,A, varargin{:});
%disp('List of all arguments:');disp(p.Results);
%disp('Using defaults for parameters:');disp(p.UsingDefaults);
r=p.Results;

maxedgewidth=7; % yEd-specific

% general parameters
nnodes = full(size(A,1));

if sum(sum(A))==0%all(all(A==0))
    fprintf('Warning: empty matrix\n');
    return;
end

% if node shape list is given => overwrite
if numel(r.nodeshapelist)
    allshapes=r.nodeshapelist;
end

% default variable names?
if numel(r.nodelabels)==0
    r.nodelabels = {};
    for i=1:nnodes
        char(64+i);
        if i<=26
            r.nodelabels{i} = char(64+i);
        else
            tmp=i-1;
            x=floor(tmp/26);
            y=mod(tmp,26)+1;
            r.nodelabels{i} = [char(64+x) char(64+y)];
        end
    end
end

% validate node colors, if given
if numel(r.nodecolors)>0
    if ~(size(r.nodecolors,1)==nnodes) || ~(size(r.nodecolors,2)==3)
        error('Node colors matrix must be %dx3', nnodes);
    end
end

% validate edge colors, if given
if numel(r.edgecolors)>0 || numel(r.edgecolorA)>0
    % both must be there
    if numel(r.edgecolors)==0
        error('If edgecolorA is given, edgecolors must also be given');
    end
    if numel(r.edgecolorA)==0
        error('If edgecolors is given, edgecolorA must also be given');
    end
    % validate size of assignments
    if ~(size(r.edgecolorA,1)==nnodes) || ~(size(r.edgecolorA,2)==nnodes)
        error('Edge color assignment matrix must be %dx%d', nnodes, nnodes);
    end
    % validate indices
    if any(r.edgecolorA(:)<0) || any(mod(r.edgecolorA(:),1)>0)
        error('Edge color assignments must be positive integer values or 0');
    end
    % validate width of color matrix
    if size(r.edgecolors,2) ~= 3
        error('Color matrix must have 3 columns');
    end
    % validate number of colors
    if max(r.edgecolorA(:)) > size(r.edgecolors,1)
        error('Edge color labels go up to %d, but color matrix is only %dx3', ...
            full(max(r.edgecolorA(:))), size(r.edgecolors,1));
    end
end

% validate groups
if numel(r.nodegroups)>0
    if ~numel(r.nodegroups)==nnodes
        error('nodegroups must be %dx1',nnodes);
    end
    ngrps=r.nodegroups;
else
    % default => all zero
    ngrps=zeros(nnodes,1);
end

% validate coordinates
if numel(r.coordinates)>0
    % must be nodes X 2
    if size(r.coordinates,2)~=2 || size(r.coordinates,1)~=nnodes
        error('Coordinates must be nodes X 2');
    end
end

% if edge weights are to be used: ensure that they are between 0 and 7
if r.useweights
    % negative values not allowed currently
    if any(A(:)<0)
        error('Negative values currently not supported');
    end
    % customscale only with scaleweights
    if numel(r.customscale)>0 && ~r.scaleweights
        error('If customscale is given, scaleweights must be true');
    end
    % scale
    if r.scaleweights
        if numel(r.customscale)>0
            % scale in given range
            minv=r.customscale(1);
            maxv=r.customscale(2);
            A=(A-minv)./(maxv-minv);
            % cut everything below and above
            A(A<0)=0;
            A(A>1)=1;
            % scale to max width
            A=A*maxedgewidth;            
        else
            % 0 and max
            A=A./max(max(A))*maxedgewidth;
        end
    end
end

% validate size of edge labels
if numel(r.edgelabels)>0
    if size(r.edgelabels,1)~=nnodes || size(r.edgelabels,2)~=nnodes
        error('edgelabels must be sparse_cell of size %dx%d', nnodes, nnodes);
    end
end

% do the shrinking, if necessary
if r.shrink
    % get involved nodes
    used=sum(A,1)~=0|(sum(A,2)~=0)';
    % do the shrinking
    A=A(used,used);
    r.nodelabels=r.nodelabels(used);
    nnodes=full(sum(used));
    if numel(r.nodecolors)>0
        r.nodecolors=r.nodecolors(used,:);
    end
    if numel(r.edgecolorA)>0
        r.edgecolorA = r.edgecolorA(used,used);
    end
    if numel(r.edgelabels)>0
        r.edgelabels=r.edgelabels(used,used);
    end
    ngrps=ngrps(used);
    if numel(r.nodeshapes)>0
        r.nodeshapes=r.nodeshapes(used);
    end
end

% open file
h=fopen(outfile,'w');
% header
writeHeader(h);
% process nodes, group-wise
ugrps=sort(unique(ngrps));

for g=1:numel(ugrps)
    % get the nodes
    curnodes=find(ngrps==ugrps(g));
    
    % if this is an actual group, prepare start
    if ugrps(g)~=0
        writeGroupHeader(h,ugrps(g));
    end
    
    % go through nodes
    for c=1:numel(curnodes)
        i=curnodes(c);
        %%% assemble node details, this is indepdentent of the groups
        % determine color
        if numel(r.nodecolors)>0
            color=colorToHTML(r.nodecolors(i,:)); % given one
        else
            color='#FFFFFF'; % white
        end
        % coords?
        if numel(r.coordinates)>0
            x=r.coordinates(i,1);
            y=r.coordinates(i,2);
        else
            x=0;
            y=0;
        end
        % determine if different shapes are used
        if numel(r.nodeshapes)>0
            shape = allshapes{r.nodeshapes(i)+1}; % +1 because it's 0-based for compatibility
        else
            shape = allshapes{1};
        end
        % write it out
        writeNode(h,i,color,r.nodelabels{i},x,y,shape);
    end
    
    
    % if this is an actual group, prepare start
    if ugrps(g)~=0
        writeGroupFooter(h,ugrps(g));
    end
    
end

% write out edges
count=0;
for e=find(A)'
    count=count+1;
    [a b] = ind2sub(size(A),e);
    if r.directed || a>b
        % define edge width/weight
        if ~r.useweights
            strwidth='1.0';
        else
            strwidth=sprintf('%.2f', full(A(a,b)));
        end
        % define edge color
        if numel(r.edgecolors)==0 || r.edgecolorA(a,b)==0
            strcolor='#000000';
        else
            strcolor=colorToHTML(r.edgecolors(r.edgecolorA(a,b),:)); % given one
        end
        % print it
        fprintf(h,'<edge id=\"e%d\" source=\"n%d\" target=\"n%d\">', count, a, b);
        fprintf(h,'<data key=\"d2\"><y:PolyLineEdge>');
        fprintf(h,'<y:LineStyle color=\"%s\" type=\"line\" width=\"%s\"/>',strcolor,strwidth);
        if ~r.directed
            fprintf(h,'<y:Arrows source=\"none\" target=\"none\"/>');
        else
            fprintf(h,'<y:Arrows source=\"none\" target=\"standard\"/>');
        end
        % label?
        if numel(r.edgelabels)>0 && numel(r.edgelabels{a,b})>0
            fprintf(h,'<y:EdgeLabel>%s</y:EdgeLabel>', r.edgelabels{a,b});
        end
        
        fprintf(h,'</y:PolyLineEdge></data></edge>\n');
    end
end
% footer
writeFooter(h);
% close file
fclose(h);

function writeNode(h,idnum,color,label,x,y,shape)
fprintf(h,'<node id= \"n%d\"><data key=\"d0\">',idnum);
fprintf(h,'<y:ShapeNode><y:Fill color=\"%s\" transparent=\"false\"/>', color);
fprintf(h,'<y:NodeLabel>%s</y:NodeLabel>', label);
if sum([x y])>0
    % also write out coordinates
    fprintf(h,'<y:Geometry x=\"%f\" y=\"%f\"/>',x,y);
end
fprintf(h,'<y:Shape type=\"%s\"/></y:ShapeNode>',shape);
fprintf(h,'</data></node>\n');

function writeGroupHeader(h,igrp)
fprintf(h,'<node id="ng%d" yfiles.foldertype="group">',igrp);
fprintf(h,'<graph id="ng%d:">',igrp);

function writeGroupFooter(h,igrp)
fprintf(h,'</graph>\n</node>');


function writeHeader(h)

fprintf(h,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(h,'<graphml xmlns="http://graphml.graphdrawing.org/xmlns/graphml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:y="http://www.yworks.com/xml/graphml" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns/graphml http://www.yworks.com/xml/schema/graphml/1.0/ygraphml.xsd">\n');
fprintf(h,'<key for="node" id="d0" yfiles.type="nodegraphics"/>\n');
fprintf(h,'<key for="edge" id="d2" yfiles.type="edgegraphics"/>\n');
fprintf(h,'<graph id="G" edgedefault="undirected">\n');

function writeFooter(h)
fprintf(h,'</graph>\n</graphml>\n');

% generate HTML color value #FFFFFF from [1.0 1.0 1.0] format
function r=colorToHTML(v)
v=floor(v*255);
r=sprintf('#%s%s%s',dec2hex(v(1),2),dec2hex(v(2),2),dec2hex(v(3),2));

function dummy
%% test on random net
labels={'A','B','C','D','E'};
grps=[2 0 1 1 2];
A = double(rand(5,5)>0.7);
A= (A+A')>0;
A=A-diag(diag(A));
writeYED('testout3.graphml',A,'nodelabels',labels,'nodegroups',grps);

%% test on sandbox multipart network
adj=g.getAdjacency;
labels=g.getLabels;
edgecolors=sparse(numel(labels),numel(labels),3);
writeYED('testout.graphml',adj.A,'nodeLabels',labels);