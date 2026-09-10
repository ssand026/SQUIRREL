function [f_result,f_string,f_inputs,didEval] = func_reduce(func, var_list, opt)
% Reduces or evaluates a function using the set of inputs "var_list"
% 
% Passes the defined variables in var_list into the inputted function and 
% attempts to evaluate the function body. If unable to evaluate, returns the 
% reduced expression for the function.
arguments
	func {mustBeA(func,{'function_handle','char','string'})}
end
arguments (Repeating)
	var_list (1,:) cell
end
arguments
	opt.varName (1,1) string = "input";
end

% check validity of var_list
if numel(var_list)==1 && iscell(var_list{1})
	var_list = var_list{1};
end

if ~all(cellfun(@numel,var_list)==2)
	error("ERROR: the input var_list has an invalid size")
elseif ~all(cellfun(@(x) isStringScalar(x{1}), var_list))
	error("ERROR: the first arguments of var_list must be a string scalar")
end

% extract the argument names/values from var_list
var_names = strings(1,numel(var_list));
var_value = cell(1,numel(var_list));
for jj = 1:numel(var_list)
	var_names(jj) = var_list{jj}{1};
	var_value{jj} = var_list{jj}{2};
end
% check for repeats
[~,indx] = unique(var_names);
var_names = var_names(indx);
var_value = var_value(indx);

% append any internalized workspace variables to the input args
if isa(func,'function_handle')
	f_work = functions(func).workspace{1};
	work_names = string(fieldnames(f_work)).';
	work_value = struct2cell(f_work).';

	% add or replace input arguments, with priority given to workspace vars
	for ii = 1:numel(work_names)
		istaken = (work_names(ii)==var_names);
		if any(istaken)
			var_value(istaken) = work_value{ii};
		else
			var_names = [var_names, work_names(ii)];
			var_value = [var_value, work_value{ii}];
		end
	end
end

% convert function to string
if isa(func,'function_handle')
	f_string = string(functions(func).function);
else
	f_string = string(func);
end
f_string = regexprep(f_string,"\s*","");

% extract possible valid variable names from the function body
f_body = regexp(f_string,"(?<=@\([^)]+\)).*","match");
[f_vars,f_notvar] = regexp(f_body,"([a-z_A-Z]+(\d+|))+","match","split");

% find variables in the function body that are defined by var_list
def_vars = intersect(f_vars,var_names);

% find variables in the function body that are not defined or built-in/existing functions
undef_vars = setdiff(f_vars,def_vars);
isafunc = arrayfun(@(var) any(exist(var)==[2,3,5]), undef_vars);
undef_vars = undef_vars(~isafunc);

% replace any defined variables with an index into the cell array f_inputs
if ~isempty(def_vars)
	input_value = cell(1,numel(def_vars));
	for ii = 1:numel(def_vars)
		f_vars(f_vars==def_vars(ii)) = opt.varName+"{"+ii+"}";
		input_value{ii} = var_value{var_names==def_vars(ii)};
	end
	% rejoin function body
	f_body = f_notvar(1) + join(reshape([f_vars; f_notvar(2:end)],1,[]),"");
end

% return the outputs
if isempty(f_body)
	% the function does not have a valid form
	f_string = "";
	f_inputs = [];
	f_result = NaN;
	didEval = false;
elseif isempty(def_vars) & isempty(undef_vars)
	% the function requires no arguments for evaluation
	f_string = "@(~) " + f_body;
	f_inputs = [];
	f_result = feval(str2func(f_string));
	didEval = true;
elseif ~isempty(def_vars) & isempty(undef_vars)
	% all of the required arguments have been defined
	f_string = "@(" + opt.varName + ") " + f_body;
	f_inputs = input_value;
	f_result = feval(str2func(f_string),input_value);
	didEval = true;
elseif ~isempty(def_vars) & ~isempty(undef_vars)
	% some of the required arguments have been defined
	f_string = "@(" + join([opt.varName,undef_vars],",") + ") " + f_body;
	f_inputs = input_value;
	f_result = [];
	didEval = false;	
elseif isempty(def_vars) & ~isempty(undef_vars)
	% none of the required arguments have been defined
	f_string = "@(" + join(undef_vars,",") + ") " + f_body;
	f_inputs = [];
	f_result = [];
	didEval = false;
end
% done
end