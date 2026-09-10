function [varargout] = meshplot3D(pts,tri,nodeData,opt)
% Visualize multiple data on a 3D mesh.
%
% If the input "nodeData" contains multiple fields (such as a matrix where each
% column is a domain-continuous coeff, or a page-array where each page is
% a simplex-continuous coeff), all data will be plotted on the same mesh and
% displayed within the same figure. It is possible to change the displayed data
% via the command fig.UserData.setIndex(...), or by using the mouse scroll wheel 
% while hovering over the figure.
%  
% SEE ALSO: PATCH, ZERO_SURFACE, MESHPLOT2D
arguments
	pts (3,:) double {mustBeFinite,mustBeReal}
	tri (:,:) {mustBeInteger,mustBePositive}
	nodeData (:,:,:) double {mustBeFinite}
	opt.plotStyle string {mustBeMember(opt.plotStyle, ...
		["surface","percent","volume","density","segment","scatter","complex"])} = "surface";
	opt.cmapLimits (1,2) {validateattributes(opt.cmapLimits,{'string','numeric'},{})} = ["min","max"];
	opt.surfValues (1,:) {validateattributes(opt.surfValues,{'string','numeric'},{})} = [];
	opt.fixedCLims (1,1) logical = false;
	opt.wrapScroll (1,1) logical = false;
	opt.titleIndex (1,1) logical = true;
	opt.axisLabels (1,1) logical = true;
	opt.edgeColor = [];
	opt.faceColor = [];
	opt.Colormap = [];
end

[Ndim, Npts] = size(pts);
tri = tri(1:Ndim+1,:);

% use the more efficient 'flat' interpolation for large meshes
if Npts > 8000
	interpMethod = "flat";
else
	interpMethod = "interp";
end

% check/prepare nodeData
%===========================================================
plotMeshOnly = false;
nodeData = squeeze(nodeData);
if isempty(nodeData)
	% only plot the mesh
	plotMeshOnly = true;
elseif ndims(nodeData)==3 | ~any(size(nodeData)==Npts)
	% invalid size for nodeData
	error("ERROR: the given node data has invalid size");
elseif size(nodeData,1)~=Npts & size(nodeData,2)==Npts
	% transpose so columns contain node-values
	nodeData = nodeData.';
end

% get the number of data frames
Nframe = size(nodeData,2);

% handle complex inputs
%-----------------------------------------------------------
% check if the complex display mode should be used
if ~isreal(nodeData) & (opt.plotStyle ~= "complex")
	warning("Data has non-real compoents. Switching to the 'complex' display-style.")
	opt.plotStyle = "complex";
end

% check if a non-complex display mode should be used
if isreal(nodeData) & (opt.plotStyle == "complex")
	warning("Data is non-complex. Switching to the 'segment' display-style.")
	opt.plotStyle = "segment";
end

if opt.plotStyle == "complex"
	opt.Colormap = hsv;
	opt.cmapLimits = (pi/2)*[-1,+1];
end


% setup colormap/limits
%===========================================================
if isnumeric(opt.cmapLimits)
	cLims = zeros(Nframe,2) + opt.cmapLimits;
end

if isnumeric(opt.surfValues) && ~isempty(opt.surfValues)
	surfVals = zeros(Nframe,size(opt.surfValues,2)) + opt.surfValues;
end

% evaluate data-dependent surface values/ cmap limits
%-----------------------------------------------------------
if isstring(opt.cmapLimits) || isstring(opt.surfValues)
	% switch between fixed and dynamic colormap limits
	if opt.fixedCLims; cData=nodeData(:); else; cData=nodeData; end

	% get statistics for nodeData
	dmin = min(cData,[],1);
	dmax = max(cData,[],1);
	dmag = max(abs(dmin),abs(dmax));
	davg = mean(cData,1);
	dstd = std(cData)+eps;

	funcEval = @(string) reshape( ...
		feval(str2func("@(min,max,mag,avg,std)"+string),dmin,dmax,dmag,davg,dstd),[],1);

	if isstring(opt.cmapLimits)
		cLims = zeros(Nframe,2);
		cLims(:,1) = funcEval(opt.cmapLimits(1));
		cLims(:,2) = funcEval(opt.cmapLimits(2));
	end

	if isstring(opt.surfValues)
		surfVals = zeros(Nframe,numel(opt.surfValues));
		for jj = 1:size(surfVals,2)
			surfVals(:,jj) = funcEval(opt.surfValues(jj));
		end
		surfVals = sort(surfVals,2);
	end
end
% ensure colormap limits are ascending (cLims(:,1) < cLims(:,2))
issame = (cLims(:,1)==cLims(:,2));
isdecr = ~(cLims(:,1)<=cLims(:,2));
cLims(issame,:) = cLims(issame,:) + eps(cLims(issame,:)).*[-1,+1];
cLims(isdecr,:) = [cLims(isdecr,2),cLims(isdecr,1)];

% set colormap if unspecified
if isempty(opt.Colormap)
	if any(sign(cLims(:,1))==-1) & any(sign(cLims(:,2))==+1)
		% use divergent colormap
		try opt.Colormap = cmap_light; end
	else
		% use monotonic colormap
		try opt.Colormap = flipud(cmap_mono); end
	end
end

% figure configs
%===========================================================
fig = figure("WindowScrollWheelFcn",@scroll_update);
fig.WindowStyle ="normal";
fig.Visible = "off";
try fig.Colormap = opt.Colormap; end

% axis settings
ax = axes("Visible","on");
ax.NextPlot = "replaceChildren";
ax.Units = "normalized";
ax.OuterPosition = [0.01 0.01 0.98 0.98];
ax.DataAspectRatio = [1 1 1];
ax.View = [-45 15];
ax.Box = "on";
ax.TickDir = "none";
try ax.Colormap = opt.Colormap; end

% axis plot-box setup
bnds = [min(pts,[],2), max(pts,[],2)];
mids = (bnds(:,1)+bnds(:,2))/2;
lims = mids + 1.05 * [-1/2, 1/2] .* (bnds(:,2)-bnds(:,1));
ax.XLim = lims(1,:);
ax.YLim = lims(2,:);
ax.ZLim = lims(3,:);

% axis labels
if opt.axisLabels
	ax.XTick = mids(1);  ax.XTickLabel = "x";
	ax.YTick = mids(2);  ax.YTickLabel = "y";
	ax.ZTick = mids(3);  ax.ZTickLabel = "z";
else
	ax.XTick = [];
	ax.YTick = [];
	ax.ZTick = [];
end

% setup initial frame
%===========================================================

% set title labeling method
%-----------------------------------------------------------
if opt.titleIndex
	% use whole indices
	label = @(x) ("   ["+sprintf("%i",x)+"]");
else
	% use fractional indices
	label = @(x) ("   ["+sprintf("%.3f",(x-1)/(Nframe-1))+"]");
end

% set the initial data index value
%-------------------------------------------------
indx = 1;
indxData.setIndex = @update_axis;
set(fig,"UserData",indxData)

% plot the mesh boundary/external faces
%-----------------------------------------------------------
[external] = tetBoundaryFacets(pts,tri);

if plotMeshOnly
	% plot boundary mesh and exit
	if isempty(opt.faceColor); opt.faceColor = "#9a8566"; end
	if isempty(opt.edgeColor); opt.edgeColor = "#000000"; end

	patch("Vertices",pts',"Faces",external,"Parent",ax, ...
		"FaceColor",opt.faceColor,"FaceAlpha",1,"EdgeColor",opt.edgeColor,"EdgeAlpha",1)
	fig.Visible = "on";
	if nargout > 0
		varargout{1} = fig;
	end
	return
else
	% plot transparent boundary mesh
	if isempty(opt.faceColor); opt.faceColor = "#242424"; end
	if isempty(opt.edgeColor); opt.edgeColor = "none"; end

	patch("Vertices",pts',"Faces",external,"Parent",ax, ...
		"FaceColor",opt.faceColor,"FaceAlpha",0.1,"EdgeColor",opt.edgeColor);
end


% interpolating contour plots
%===========================================================
if any(opt.plotStyle==["surface","percent","volume"])
	spts = cell(1,Nframe);
	surfs = cell(1,Nframe);
	colorData = cell(1,Nframe);

	if opt.plotStyle=="percent" && opt.fixedCLims
		% use fixed surface values
		shells = [0.023, 0.067, 0.159, 0.309, 0.500]; % cummulative z-scores
		[vals,~] = enclosingSurface(nodeData(:),shells,0,"outerSplit");
		[fixvals,~] = unique(vals);
	end


	% get interpolated surfaces
	for jj = 1:Nframe
		data = nodeData(:,jj);

		if all(abs(data-data(1)) < 10*eps(data(1)))
			% all data is same-valued
			continue
		elseif opt.plotStyle == "surface"
			% contour plot: evenly-spaced shells OR custom surface values
			%-----------------------------------------------------------
			Nshells = 6;
			if isempty(opt.surfValues)
				limSign = sign(cLims(jj,:));
				if isequal(limSign,[-1,+1])
					% log spacing
					r = ceil(Nshells*abs(cLims(jj,:))./diff(cLims(jj,:)));
					log2x = @(a) (linspace((1-sign(a))/2,(1+sign(a))/2,abs(a)+1));
					vals = [cLims(jj,1).*log2x(r(1)), cLims(jj,2).*log2x(r(2))];
				elseif isequal(limSign,[-1,0]) | isequal(limSign,[0,+1])
					% equal spacing plus zero point
					vals = linspace(cLims(jj,1),cLims(jj,2),Nshells+1);
				elseif isequal(limSign,fliplr(limSign))
					% equal spacing
					vals = linspace(cLims(jj,1),cLims(jj,2),Nshells);
				else
					continue
				end
				vals(vals==0) = [];
			else
				vals = surfVals(jj,:);
			end

		elseif opt.plotStyle == "percent"
			% shells encapsulate fixed percentage of data
			%-----------------------------------------------------------
			if opt.fixedCLims
				% use fixed data shells
				vals = fixvals;
			else
				% get portion of data encapsulated by each contour shell
				shells = [0.023, 0.067, 0.159, 0.309, 0.500]; % cummulative z-scores
				[vals,~] = enclosingSurface(data,shells,0,"outerSplit");
				vals = unique(vals);
			end
			% omit surfaces that are too close to zero
			dcut = 0.075;
			islarge = (abs(vals) > dcut*max(abs(data)));
			vals = vals(islarge);

		elseif opt.plotStyle == "volume"
			% shells encapsulate fixed volumes
			%-----------------------------------------------------------
			shells = [0.01, 0.025, 0.05 : 0.15 : 0.95, 0.975, 0.99].';
			tavg = mean(data(tri),1);
			[tavg,iavg] = sort(tavg);
			svol = simplex_vol(pts,tri);
			tvol = cumsum(svol(iavg));
			tvol = tvol/tvol(end);

			vals = zeros(length(shells),1);
			for ll = 1:length(shells)
				[a] = find(tvol < shells(ll),1,"last");
				[b] = find(tvol > shells(ll),1,"first");
				try 
					vals(ll) = abs([tvol(a),tvol(b)]-shells(ll)) / ...
					abs(tvol(b)-tvol(a))*[tavg(b);tavg(a)];
				catch
					vals(ll) = 0;
				end
			end
			% omit surfaces that are too close to zero
			dcut = 0.05;
			islarge = (abs(vals) > dcut*max(abs(vals)));
			vals = vals(islarge);
		end

		% remove any simplices with all points outside the surfval limits
		pos = vals(vals>=0); if isempty(pos); pos = 0; end
		neg = vals(vals<=0); if isempty(neg); neg = 0; end
		inner = (data > max(neg)) & (data < min(pos));
		lower = (data < min(neg)); upper = (data > max(pos));
		isOutside = all(inner(tri),1) | all(lower(tri),1) | all(upper(tri),1);
		tSubset = tri(:,~isOutside);

		% find the interpolated surface for each surfval
		if ~isempty(vals)
			for kk = 1:length(vals)
				[newpts,newfaces] = zero_surface(pts,tSubset,data - vals(kk));
				surfs{jj} = [surfs{jj}; newfaces+size(spts{jj},1)];
				spts{jj} = [spts{jj}; newpts];
				colorData{jj} = [colorData{jj}; kk * ones(size(newfaces,1),1)];
			end
			colorData{jj} = vals(colorData{jj});
		end
	end

	% plot surfaces
	plt = patch("Vertices",spts{indx},"Faces",surfs{indx},"Parent",ax, ...
		"FaceColor","flat","EdgeColor","none","FaceAlpha",0.3);

	% set colordata
	set(plt,"Cdata",colorData{indx});
end

% density/transparent style plot
%===========================================================
if opt.plotStyle == "density"
	[allFaces] = mesh_faces(pts,tri);
	[internal] = setdiff(allFaces,external,"rows");

	% plot internal faces
	plt = patch("Vertices",pts.',"Faces",internal,"Parent",ax, ...
		"FaceColor","interp","FaceAlpha",interpMethod,"EdgeColor","none");

	% set color map
	colorData = nodeData;
	set(plt,"Cdata",colorData(:,indx));

	% set alpha map
	alphaData = abs(nodeData);
	set(plt,"AlphaDataMapping","scaled")
	set(plt,"FaceVertexAlphaData",alphaData(:,indx))
end

% segment/wireframe style plot
%===========================================================
if opt.plotStyle == "segment"
	[allFaces] = mesh_faces(pts,tri);
	[internal] = setdiff(allFaces,external,"rows");

	% plot mesh edges
	plt=patch("Vertices",pts.',"Faces",internal,"Parent",ax, ...
		"FaceColor","none","FaceAlpha",0,"EdgeColor","interp","EdgeAlpha","flat","lineWidth",3);

	% set color map
	colorData = nodeData;
	set(plt,"Cdata",colorData(:,indx));

	% set alpha map
	alphaData = abs(nodeData);
	set(plt,"AlphaDataMapping","scaled");
	set(plt,"FaceVertexAlphaData",alphaData(:,indx));
end

% scatter/points style plot
%===========================================================
if opt.plotStyle == "scatter"

	% permit display of scatter plot
	ax.NextPlot = "add";

	% get marker size data
	sizeData = abs(nodeData);
	sbnds = mean(sizeData,1) + std(sizeData,0,1).*[-1;2]; sbnds(sbnds<0)=0;
	sizeData = (sizeData - sbnds(1,:)) ./ (sbnds(2,:)-sbnds(1,:));
	sizeData(sizeData<0) = 0;
	sizeData = 1 + 24*(sizeData);

	% set color map
	colorData = nodeData;

	% plot points/dots
	[x,y,z] = deal(pts(1,:).',pts(2,:).',pts(3,:).');
	plt=scatter3(x,y,z,sizeData(:,indx),colorData(:,indx),"Parent",ax, ...
		"Marker","square","MarkerFaceColor","flat","MarkerFaceAlpha",1,"LineWidth",1);

	ax.NextPlot = "replacechildren";
end

% phase-space style plot
%===========================================================
if opt.plotStyle == "complex"
	[allFaces] = mesh_faces(pts,tri);
	[internal] = setdiff(allFaces,external,"rows");

	% plot mesh edges
	plt=patch("Vertices",pts.',"Faces",internal,"Parent",ax, ...
		"FaceColor","none","FaceAlpha",0,"EdgeColor",interpMethod,"EdgeAlpha","flat","lineWidth",3);

	% set color map
	colorData = angle(nodeData); % phase information
	set(plt,"Cdata",colorData(:,indx));

	% set alpha map
	alphaData = abs(nodeData).^2; % magnitude information
	set(plt,"AlphaDataMapping","scaled");
	set(plt,"FaceVertexAlphaData",alphaData(:,indx));
end


% display figure and define update functions
%===========================================================
set(fig,"Name",label(indx)); % set figure title
set(ax,"CLim",cLims(indx,:)); % set figure clims
fig.Visible = "on";
if nargout > 0
	varargout{1} = fig;
end

%-----------------------------------------------------------
	function scroll_update(~,event)
	% use mouse scroll wheel to update the displayed data index
	new_indx = indx + (event.VerticalScrollCount);
	update_axis(new_indx)
	end
%-----------------------------------------------------------
	function update_axis(new_indx)
	% updates the displayed data after performing checks
	doUpdate = checkIndex(new_indx);

	if ~doUpdate; return; end
	% update style-specific data
	switch opt.plotStyle
		case {"surface","percent","volume"}
			set(plt,"CData",colorData{indx});
			set(plt,"Vertices",spts{indx});
			set(plt,"Faces",surfs{indx});
		case "scatter"
			set(plt,"SizeData",sizeData(:,indx));
			set(plt,"CData",colorData(:,indx));
		otherwise
			set(plt,"CData",colorData(:,indx))
			set(plt,"FaceVertexAlphaData",alphaData(:,indx));
	end
	set(fig,"Name",label(indx)); % update figure title
	set(ax,"CLim",cLims(indx,:)); % update color-data limits
	end
%-----------------------------------------------------------
	function [doUpdate] = checkIndex(new_indx)
	% validates new_indx, enforces indexing bounds, and determines if update is needed
	doUpdate = false;
	if isnumeric(new_indx) && (round(new_indx)==new_indx)
		if opt.wrapScroll
			new_indx = mod(new_indx, Nframe);
			if new_indx==0; new_indx = Nframe; end
		else
			new_indx = min(new_indx,Nframe);
			new_indx = max(new_indx,1);
		end
		if new_indx~=indx
			indx = new_indx;
			doUpdate = true;
		end
	end
	end
%-----------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dvals,dfrac] = enclosingSurface(data,dataFrac,midpoint,surfType)
% finds the surface that encloses a given fraction of the dataset.
arguments
	data   (:,1)
	dataFrac (:,1) {mustBeInRange(dataFrac,0,1)}
	midpoint (1,1) double = 0;
	surfType string {mustBeMember(surfType, ...
		["inner","outer","innerSplit","outerSplit","absolute"])} = "absolute";
end
% misc
Nbin = 256;
scl = max(abs(data));
tol = 4*eps(scl);

% bin data and find the relative frequency of the data values
binEdge = [linspace(tol, scl, Nbin/2-1), scl+tol];
binEdge = [-fliplr(binEdge), 0, binEdge];
bins = (binEdge(1:end-1)+binEdge(2:end))/2;
freq = histcounts(data,binEdge);
bins = bins(:); freq = freq(:);

% find the data to the left and right of the data midpoint
[~,cntr_indx] = mink(abs(bins-midpoint),2);
indL = [1:min(cntr_indx)];
indR = [max(cntr_indx):Nbin];

% if surfType ~= "absolute"
% 	freq(abs(bins) < tol) = 0; % ignore zeros
% end
freq = freq / sum(freq);

switch surfType
	case {"inner", "outer"}
		% fold data values across the data midpoint
		bins = [midpoint-bins(indL); bins(indR)-midpoint];
		[bins,sort_indx] = sort(bins);
		freq = freq(sort_indx);

		% combine overlapping bins
		[~,unqL] = unique(bins,"first");
		[~,unqR] = unique(bins,"last");
		bins = bins(unqL);
		freq = freq(unqL) + (unqL~=unqR).*freq(unqR);

	case {"innerSplit", "outerSplit"}
		% normalize the data distribution for both sides
		lsum = sum(freq(indL));
		rsum = sum(freq(indR));
		lfreq = freq(indL)/(lsum + ~lsum);
		rfreq = freq(indR)/(rsum + ~rsum);
end

switch surfType
	case "absolute"
		freq = cumsum(freq);
		[dvals,dfrac] = enclosingVals(bins,freq,dataFrac,"lower");
	case "inner"
		[dvals,dfrac] = enclosingVals(bins,freq,dataFrac,"lower");
		dvals = midpoint + [-dvals, +dvals];
		dfrac = [-dfrac, dfrac];
	case "outer"
		[dvals,dfrac] = enclosingVals(bins,freq,dataFrac,"upper");
		dvals = midpoint + [-dvals, +dvals];
		dfrac = [-dfrac, dfrac];
	case "innerSplit"
		[lvals,lfrac] = enclosingVals(bins(indL),lfreq,dataFrac,"upper");
		[rvals,rfrac] = enclosingVals(bins(indR),rfreq,dataFrac,"lower");
		dvals = [lvals, rvals];
		dfrac = [-lfrac, rfrac];
	case "outerSplit"
		[lvals,lfrac] = enclosingVals(bins(indL),lfreq,dataFrac,"lower");
		[rvals,rfrac] = enclosingVals(bins(indR),rfreq,dataFrac,"upper");
		dvals = [lvals, rvals];
		dfrac = [-lfrac, rfrac];
end
% final sort
[dvals,sort_indx] = sort(dvals);
dfrac = dfrac(sort_indx);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataval,dataFrac] = enclosingVals(data,freq,dataFrac,dataPart)
% finds data values that enclose the given fraction "dataFrac" of
% the dataset defined by a list of values "data", that appear in
% the dataset with relative frequency "freq". The input "dataPart"
% determines if dataval encloses the lower or upper part of the data
arguments
	data (1,:) double
	freq (1,:) double
	dataFrac (1,:) double {mustBeInRange(dataFrac,0,1)}
	dataPart {mustBeMember(dataPart,["lower","upper"])}
end

% sort the data
[data,sort_indx] = sort(data);
freq = freq(sort_indx);
tol = eps(max(abs(data)));

% find the cummulative data value distributions
if dataPart=="lower"; freq = cumsum(freq,"forward"); end
if dataPart=="upper";  freq = cumsum(freq,"reverse"); end

% eliminate repeated entries from the frequency distribution
[~,indL] = unique(freq,"first");
[~,indR] = unique(freq,"last");
data = [data(indL)-2*tol, data(indR)+2*tol];
freq = [freq(indL)-2*eps, freq(indR)+2*eps];

[~,sort_indx] = sort(data);
data = data(sort_indx);
freq = freq(sort_indx);

% find the surface values that enclose the specified proportions of the data
dataFrac = unique(dataFrac);
dataval = interp1(freq,data,dataFrac);

% only return valid surfaces
isvalid = ~isnan(dataval);
dataval = dataval(isvalid);
dataFrac = dataFrac(isvalid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
