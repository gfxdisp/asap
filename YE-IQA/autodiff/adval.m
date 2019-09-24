function y = adval(df,varargin)

%ADVAL Internal processing of derivative data.

% df       - function handle of derivative
% varargin - {list of arguments, index of variable}


% default if only one argument
if numel(varargin)==1
  varargin{2} = 1;
end

% prepare variable and result
ivar = varargin{end};    % index of variable
x = varargin{ivar};      % variable
arg = varargin(1:end-1); % all arguments
y = x;                   % allocate result

% exit if x empty
if isempty(x)
  return
end

% determine calling function fct
dbs = dbstack;
[~,fct] = dbs.name;
l = find(fct=='.',1,'last');
if ~isempty(l)
  fct = fct(l+1:end);
end

% convert audi array to singleton
if numel(x)>1
  sc = size(x(1).c);                   % size of data stored in .c
  z = x(1);                            % tmp audi singleton
  z.c = zeros(sc(1)*numel(x),sc(2));   % allocate .c array
  I = 1:sc(1);
  for i = 1:numel(x)                   % vertical catenation of array elements
    z.c(I,:) = x(i).c;
    I = I + sc(1);
  end
  x = z;                               % overwrite x
end

% compute audi singleton f = fct(x)
f = x;      % initialize result
if x.k==0 || (x.n==1 && all(all(x.c(:,2:end)==[1,zeros(1,x.k-1)])))
  % values only (k==0) or native Taylor expansion of fct(x)
  x0 = x.c(:,1);         % point of evaluation
  arg{ivar} = x0;        % convert list of arguments to double
  % evaluate f0 = fct(x0)
  try
    f0 = feval(fct,arg{:});
  catch
    f0 = feval(['audi.' fct],arg{:});
  end
  if x.k==0              % values only
    f.c = f0;
  else                   % native Taylor
    xx= x;               % xx = variable with reduced order of differentiation
    xx.k = xx.k - 1;
    xx.c = xx.c(:,1:xx.k+1);
    arg{ivar} = xx;
    % derivative dff = df(xx)
    try
      dff = feval(df,arg{:});
    catch
      dff = feval(['audi.' df],arg{:});
    end
    f.c  = [f0 dff.c];   % catenate values and derivatives
  end
else
  % general case: determine native Taylor polyomial and evaluate at x
  t = ainit(x.c(:,1),x.k);     % variable for native Taylor polynomial a
  t.s = x.s;
  arg{ivar} = t;
  x.c(:,1) = 0;                % remainder, f = f0 + a(x)
  % native Taylor polynomial sum_i a(i)/i! (t-t0)^i
  try
    a = feval(fct,arg{:});     
  catch
    a = feval(['audi.' fct],arg{:});
  end
  a.c = a.c./[1 cumprod(1:a.k)];
  f = 0.*x;
  f.c(:,1) = a.c(:,end);
  for i = a.k:-1:1    % evaluation by Horner
    f = x.*f;
    f.c(:,1) = a.c(:,i);
  end
end

% reshape singleton f to obtain array y
if numel(y)>1
  I = 1:sc(1);
  for i = 1:numel(y)
    y(i).c = f.c(I,:);
    I = I + sc(1);
  end
else
  y = f;
end