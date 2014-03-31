% TODO: write docu
function writeCytoscape(outfile,varargin)
% TODO: edge labels buggy

% pre-defined shapes
allshapes={'ellipse','diamond','rectangle','triangle'};

maxedgewidth=12; % cytoscape-specific

% if the first argument is not a structure, it's A
if ~isstruct(varargin{1})
    A=varargin{1};
    varargin=varargin(2:end);
else
    A=varargin{1}.A;
end

% parse parameters
p = inputParser;
p.StructExpand=true; % if we allow parameter  structure expanding
p.addRequired('outfile', @ischar); % ok
p.addRequired('A', @isnumeric); % ok
p.addParamValue('nodelabels', '', @iscellstr); % ok
p.addParamValue('nodecolors',[],@isnumeric); % ok
p.addParamValue('nodeannos',[],@isstruct); % ok
% p.addParamValue('nodegroups',[],@isnumeric); % not now
p.addParamValue('directed', 0, @islogical); % ok
p.addParamValue('shrink', true, @islogical); % ok
p.addParamValue('useweights', false, @islogical); % ok
p.addParamValue('scaleweights', true, @islogical); % ok
p.addParamValue('customscale', [], @isnumeric); % ok
p.addParamValue('edgecolors',[],@isnumeric); % ok
p.addParamValue('edgecolorA',[],@isnumeric); % ok
% p.addParamValue('coordinates',[],@isnumeric); % not now
p.addParamValue('edgelabels',[],@(x)isa(x,'sparse_cell')); % currently buggy
p.addParamValue('edgeannos',[],@isstruct); % ok
p.addParamValue('nodeshapes',[],@isnumeric);
p.addParamValue('nodeshapelist',[],@iscellstr);
p.addParamValue('nodesizes',[],@isstruct);
p.addParamValue('nodeemptylabel',[],@islogical);
p.parse(outfile,A, varargin{:});
%disp('List of all arguments:');disp(p.Results);
%disp('Using defaults for parameters:');disp(p.UsingDefaults);
r=p.Results;

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
    error('edgelabels is currently buggy for cytoscape');
    if size(r.edgelabels,1)~=nnodes || size(r.edgelabels,2)~=nnodes
        error('edgelabels must be sparse_cell of size %dx%d', nnodes, nnodes);
    end
end

% validate node size information
if numel(r.nodesizes)>0
    % validate structure
    if ~isfield(r.nodesizes,'width') || ~isfield(r.nodesizes,'height')
        error('nodesizes must be a structure with ''width'' and ''height'' fields');
    end
    % validate field dimensions
    if numel(r.nodesizes.width)~=nnodes
        error('width field must contain %d values', nnodes);
    end
    if numel(r.nodesizes.height)~=nnodes
        error('height field must contain %d values', nnodes);
    end
end

% validate empty label info
if numel(r.nodeemptylabel)>0
    if numel(r.nodeemptylabel)~=nnodes
        error('nodeemptylabel must contain %d logical values', nnodes);
    end
end

% validate node annotations
nodeannoTypes=[];
if numel(r.nodeannos)>0
    % verify all have a supported type and have the correct length
    f=fields(r.nodeannos);
    for i=1:numel(f)
        % check type
        if ~iscell(r.nodeannos.(f{i})) && ~isnumeric(r.nodeannos.(f{i}))
            error('Unsupported type for annotation field ''%s''', f{i});
        end
        % save specific type information
        if iscell(r.nodeannos.(f{i}))
            nodeannoTypes.(f{i})='string';
        elseif isnumeric(r.nodeannos.(f{i}))
            % integers?
            if all(mod(r.nodeannos.(f{i})(~isnan(r.nodeannos.(f{i}))),1)==0)
                % integer
                nodeannoTypes.(f{i})='integer';
            else
                % some floating point numbers
                nodeannoTypes.(f{i})='real';
            end
        end
        % check length
        if numel(r.nodeannos.(f{i})) ~= nnodes
            error('node annotation field ''%s'' must contain %d elements', f{i}, nnodes);
        end
    end
end

% validate edge annotations
edgeannoTypes=[];
if numel(r.edgeannos)>0
    % verify all have a supported type and have the correct length
    f=fields(r.edgeannos);
    for i=1:numel(f)
        % check type
        if ~iscell(r.edgeannos.(f{i})) && ~isnumeric(r.edgeannos.(f{i}))
            error('Unsupported type for annotation field ''%s''', f{i});
        end
        % save specific type information
        if iscell(r.edgeannos.(f{i}))
            edgeannoTypes.(f{i})='string';
        elseif isnumeric(r.edgeannos.(f{i}))
            % integers?
            if all(all(mod(r.edgeannos.(f{i})(~isnan(r.edgeannos.(f{i}))),1)==0))
                % integer
                edgeannoTypes.(f{i})='integer';
            else
                % some floating point numbers
                edgeannoTypes.(f{i})='real';
            end
        end
        % check length
        if size(r.edgeannos.(f{i}),1) ~= nnodes || size(r.edgeannos.(f{i}),2) ~= nnodes
            error('edge annotation field ''%s'' must be contain %d X %d elements', f{i}, nnodes, nnodes);
        end
    end
end



% do the shrinking, if necessary
if r.shrink
    % get involved nodes
    used=sum(A,1)>0|(sum(A,2)>0)';
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
    % also shrink annotations
    if numel(r.nodeannos)
        f=fields(r.nodeannos);
        for i=1:numel(f)
            r.nodeannos.(f{i})=r.nodeannos.(f{i})(used);
        end
    end
    if numel(r.edgeannos)
        f=fields(r.edgeannos);
        for i=1:numel(f)
            r.edgeannos.(f{i})=r.edgeannos.(f{i})(used,used);
        end
    end
    if numel(r.nodeshapes)>0
        r.nodeshapes=r.nodeshapes(used);
    end
    if numel(r.nodesizes)>0
        r.nodesizes.width=r.nodesizes.width(used);
        r.nodesizes.height=r.nodesizes.height(used);
    end
    if numel(r.nodeemptylabel)>0
        r.nodeemptylabel=r.nodeemptylabel(used);
    end
end



% open file
h=fopen(outfile,'w');
% header
fprintf(h,'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n');
fprintf(h,'<graph label="Network" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML"  directed="%d">\n', r.directed);
fprintf(h,'<att name="documentVersion" value="1.1"/>\n');


%%%% NODES
for i=1:size(A,1)
    % determine color
    if numel(r.nodecolors)>0
        color=colorToHTML(r.nodecolors(i,:)); % given one
    else
        color='#FFFFFF'; % white
    end
    % determine if different shapes are used
    if numel(r.nodeshapes)>0
        shape = allshapes{r.nodeshapes(i)+1}; % +1 because it's 0-based for compatibility
    else
        shape = allshapes{1};
    end
    % width and height
    sizestr='';
    if numel(r.nodesizes)>0
        if r.nodesizes.width(i)>0
            sizestr = sprintf('%sw=''%.2f'' ', sizestr, r.nodesizes.width(i));
        end
        if r.nodesizes.height(i)>0
            sizestr = sprintf('%sh=''%.2f'' ', sizestr, r.nodesizes.height(i));
        end
    end       
    % write it out
    label=r.nodelabels{i};
    fprintf(h,'<node label="%s" id="%d">\n',label,i);
    fprintf(h,'<att type="string" name="NODE_TYPE" value="DefaultNode"/>\n');
    fprintf(h,'<graphics type="%s" cy:nodeLabel="%s" fill="%s" %s/>\n',shape,label,color,sizestr);
    % empty label?
    if numel(r.nodeemptylabel)>0 && r.nodeemptylabel(i)==true
        fprintf(h, '<att type="string" name="node.label" value=""/>\n');
    end
    % any further annotation?
    if numel(r.nodeannos)>0
        f=fields(r.nodeannos);
        for j=1:numel(f)
            toprint='';
            if strcmp(nodeannoTypes.(f{j}),'string')
                % cell entry
                cur=r.nodeannos.(f{j}){i};
                if numel(cur)>0
                    if ischar(cur)
                        toprint=cur;
                    else
                        % just skip
                      %  error('There is an unsupported entry type in ''%s''',f{j});
                    end
                end
            elseif strcmp(nodeannoTypes.(f{j}),'real')
                % numeric
                if ~isnan(r.nodeannos.(f{j})(i))
                    toprint=sprintf('%f',r.nodeannos.(f{j})(i));
                end
            elseif strcmp(nodeannoTypes.(f{j}),'integer')
                % numeric
                if ~isnan(r.nodeannos.(f{j})(i))
                    toprint=sprintf('%d',r.nodeannos.(f{j})(i));
                end
            else
                error('Bug! Argument control does not work');
            end
            % if we have something to print, do it now
            if numel(toprint)>0
                fprintf(h,'<att type="%s" name="%s" value="%s"/>\n', nodeannoTypes.(f{j}), f{j}, XMLEncode(toprint));
            end
        end
    end
    
    
    fprintf(h,'</node>\n');
    
    
end


%%%% EDGES

% write out edges
count=0;
for e=find(A)'
    count=count+1;
    [a b] = ind2sub(size(A),e);
    if r.directed || a>b
        % define edge width/weight
        if ~r.useweights
            strwidth='3.0';
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
        fprintf(h,'<edge label=\"%s (EDGEMARKER) %s\" source=\"%d\" target=\"%d\">\n', r.nodelabels{a},r.nodelabels{b}, a, b);
        %         fprintf(h,'<data key=\"d2\"><y:PolyLineEdge>');
        %         fprintf(h,'<y:LineStyle color=\"%s\" type=\"line\" width=\"%s\"/>',strcolor,strwidth);
        if ~r.directed
            targetArrow='0';
        else
            targetArrow='1';
        end
        % label?
        if numel(r.edgelabels)>0 && numel(r.edgelabels{a,b})>0
            labelstr=sprintf('cy:edgeLabel="%s" ', r.edgelabels{a,b});
        else
            labelstr='';
        end
        % write it out
        fprintf(h,'<graphics %swidth="%s" fill="%s" cy:sourceArrow="0" cy:targetArrow="%d" />\n', labelstr, strwidth, strcolor, targetArrow);
        
        % any further annotation? (copied from node anno writer above
        if numel(r.edgeannos)>0
            f=fields(r.edgeannos);
            for j=1:numel(f)
                toprint='';
                if strcmp(edgeannoTypes.(f{j}),'string')
                    % cell entry
                    cur=r.edgeannos.(f{j}){a,b};
                    if numel(cur)>0
                        if ischar(cur)
                            toprint=cur;
                        else
                            % just skip
                          %  error('There is an unsupported entry type in ''%s''',f{j});
                        end
                    end
                elseif strcmp(edgeannoTypes.(f{j}),'real')
                    % numeric
                    if ~isnan(r.edgeannos.(f{j})(a,b))
                        toprint=sprintf('%f',r.edgeannos.(f{j})(a,b));
                    end
                elseif strcmp(edgeannoTypes.(f{j}),'integer')
                    % numeric
                    if ~isnan(r.edgeannos.(f{j})(a,b))
                        toprint=sprintf('%d',r.edgeannos.(f{j})(a,b));
                    end
                else
                    error('Bug! Argument control does not work');
                end
                
                % if we have something to print, do it now
                if numel(toprint)>0
                    fprintf(h,'<att type="%s" name="%s" value="%s"/>\n', edgeannoTypes.(f{j}), f{j}, XMLEncode(toprint));
                end
            end
        end
        
        
        fprintf(h,'</edge>\n');
        
    end
end



% footer
fprintf(h,'</graph>');
% close file
fclose(h);



% generate HTML color value #FFFFFF from [1.0 1.0 1.0] format
function r=colorToHTML(v)
v=floor(v*255);
r=sprintf('#%s%s%s',dec2hex(v(1),2),dec2hex(v(2),2),dec2hex(v(3),2));


function dummy
%% test code
clc;
labels={'XA','B','AAC','D','E2'};
grps=[2 0 1 1 2];
A = zeros(5);
A(1,2)=1;
A(1,3)=2;
A(2,4)=3;
A(3,4)=4;
A(4,5)=2;
edgelabels=sparse_cell(5,5);
% edgelabels{4,2}='a label';
% edgelabels{2,1}='p=2x10e-4';
A=A+A';
A=A-diag(diag(A));
edgeColors=[0 0 0 ; 1 0 0];
edgeColorA=zeros(5);
edgeColorA(A~=0)=randi(2,sum(A(:)~=0),1);

% construct some node annotations
nodeannos=[];
nodeannos.type=cell(5,1);
nodeannos.type{2}='Amino acid';
nodeannos.type{4}='Amino acid';
nodeannos.type{5}='Lipid';
nodeannos.pval=rand(5,1);

edgeannos.pval=rand(5,5);

writeCytoscape('test.xgmml',A, 'nodelabels',labels, 'nodecolors', rand(5,3), 'directed', false, ...
    'edgecolors', edgeColors, 'edgecolora', edgeColorA, 'useweights', true, 'scaleweights', true,...
    'nodeannos', nodeannos, 'edgeannos', edgeannos );
fprintf('Done.\n');