classdef audi
  %AUDI MATLAB class for automatic differentiation.
  %   The AUDI class establishes automatic differentiation by means of
  %   operator overloading. It provides
  %
  %   * Ordinary and partial derivatives of arbitrary order
  %   * Differential operators, like Laplacian and curl
  %   * Taylor expansions of explicit expressions and solutions of ODEs
  %   * Curvature computation for curves and surfaces
  %   * Resolution of singularities of type "0/0" by l'Hospital's rule
  %
  %   Detailed information can be found in the <a href="matlab: ahelp(0)">Documentation</a>
  %   and the <a href="matlab: ahelp(-1)">Examples</a>.
  %
  %   Version: 3.1
  %   Date:    March 31, 2017
  %   Author:  Ulrich Reif
  
  %% PROPERTIES
  properties
    n       % number of variables
    k       % maxial order of differentiation
    c       % matrix of derivative data
    s       % size of data
  end
  %% METHODS
  methods
    %% CONSTRUCTOR
    function r = audi(a,d,n,k)
      if nargin==0                             % empty audi
        r.c = [];
        return
      end
      if k<0
        error('Order must not be negative.')
      end
      audi.mat(k,n);
      z = zeros(1,n);                          % index vector of zeros
      r.n = n;
      r.k = k;
      r.s = size(a);
      r.c = zeros(prod(r.s),nchoosek(n+k,k));
      r.c(:,1) = a(:);
      if d>0 && k>0
        z(d) = 1;
        r.c(:,idx(r,z)) = 1;
      end
      sumidx(r);
    end
    %% USER-DEFINED OVERLOADS
    % Functions you want to overload can be entered here.
    %
    % function y = fct(x)
    %   df = @(x) ...
    %   y = adval(df,x);
    % end
    %% ELEMENTARY FUNCTIONS
    function y = sin(x)
      df = @(x) cos(x);
      y = adval(df,x);
    end
    function y = asin(x)
      df = @(x) (1 - x.^2).^(-1/2);
      y = adval(df,x);
    end
    function y = cos(x)
      df = @(x) -sin(x);
      y  = adval(df,x);
    end
    function y = acos(x)
      y = pi/2 - asin(x);
    end
    function y = exp(x)
      df = @(x) exp(x);
      y  = adval(df,x);
    end
    function y = log(x)
      df = @(x) 1./x;
      y  = adval(df,x);
    end
    function y = tan(x)
      df = @(x) 1 + tan(x).^2;
      y  = adval(df,x);
    end
    function y = atan(x)
      df = @(x) 1./(1 + x.^2);
      y  = adval(df,x);
    end
    function y = atan2(x2,x1)
      y = angle(x1 + 1i*x2);
    end
    function y = sinh(x)
      df = @(x) cosh(x);
      y  = adval(df,x);
    end
    function y = asinh(x)
      df = @(x) (x.^2 + 1).^(-1/2);
      y  = adval(df,x);
    end
    function y = cosh(x)
      df = @(x) sinh(x);
      y  = adval(df,x);
    end
    function y = acosh(x)
      df = @(x) (x.^2 - 1).^(-1/2);
      y  = adval(df,x);
    end
    function y = tanh(x)
      df = @(x) 1 - tanh(x).^2;
      y  = adval(df,x);
    end
    function y = atanh(x)
      df = @(x) 1./(1 - x.^2);
      y  = adval(df,x);
    end
    function y = sqrt(x)
      y = x.^(1/2);
    end
    function y = abs(x)
      y = x;
      for i = 1:numel(x)
        if isreal(x(i))
          y(i) = -min(x(i),-x(i));
        else
          y(i) = sqrt(x(i).*conj(x(i)));
        end
      end
    end
    function y = sign(x)
      df = @(x) 0*x;
      y  = adval(df,x);
    end
    %% POLYNOMIALS
    function y = polyval(a,x)
      y = 0.*x + a(1);
      for i = 2:numel(a)
        y = y.*x + a(i);
      end
    end
    function y = polyvalm(a,x)
      e = eye(size(x));
      y = 0.*x + a(1)*e;
      for i = 2:numel(a)
        y = y*x + a(i)*e;
      end
    end
    function y = roots(a)
      if ~isvector(a)
        error('input must be vector')
      end
      na = numel(a);
      if a(1)==0
        y = roots(a(2:na));
        return
      elseif a(na)==0
        y = [0.*a(1); roots(a(1:na-1))];
        return
      else
        a = a./a(1);
      end
      if na <= 1
        y = [];
      elseif na==2
        y = -a(2);
      elseif na==3
        p = (-.5).*a(2);
        d = sqrt(p.*p-a(3));
        y = [p-d;p+d];
      elseif na==4
        p = a(3) - (1/3).*a(2).*a(2);
        q = a(2).*((2/27).*a(2).*a(2)-(1/3).*a(3))+a(4);
        w = -.5+.5i*sqrt(3);
        t = ((-.5).*q+sqrt(.25.*q.*q+(1/27).*p.*p.*p)).^(1/3)*[1;w;w'];
        y = t - (1/3).*(p./t+a(2));
        if isreal(a)
          for i = 1:3
            j = abs(imag(y(i).c)./real(y(i).c))<1e-12;
            y(i).c(j) = real(y(i).c(j));
          end
        end
      elseif size(a(1).c,1)==1
        a0 = aeval(a,0);
        r0 = roots(a0);
        g0 = 1./polyval(polyder(a0),r0);
        y = r0;
        for i=1:a(1).k
          y = y - polyval(a,y).*g0;
        end
        y = y - aeval(y,0) + r0;
      else
        error('for degree greater than 3, coefficients must be scalar')
      end
    end
    function y = polyder(a)
      na = size(a,2) - 1;
      y = (na:-1:1).*a(:,1:na);
    end
    function y = polyint(a,k)
      if nargin==1
        k = 0;
      end
      na = size(a,2);
      y = [(na:-1:1).\a(:,1:na) k+0.*a(:,1)];
    end
    function a = polyfit(x,y,n)
      A(:,n+1) = 0.*x(:) + 1;
      for i = n:-1:1
        A(:,i) = x(:).*A(:,i+1);
      end
      a = (A\y(:)).';
    end
    function y = conv(a,b)
      if numel(a)<numel(b)
        y = conv(b,a);
      else
        y = repmat(0.*a(1).*b(1),numel(a)+numel(b)-1,1);
        j = 0:numel(a)-1;
        for i = 1:numel(b)
          y(i+j) = y(i+j) + a(:).*b(i);
        end
        if isrow(a)
          y = y.';
        end
      end
    end
    function [y,r] = deconv(a,b)
      r = a(:);
      b = b(:);
      ny = numel(a)-numel(b)+1;
      y = repmat(0.*a(1).*b(1),ny,1);
      j = 0:numel(b)-1;
      for i = 1:ny
        y(i) = r(i)./b(1);
        r(i+j) = r(i+j) - y(i).*b;
      end
      if isrow(a)
        r = r.';
        y = y.';
      end
    end
    %% VECTORS AND MATRICES
    function y = cross(a,b)
      y = a;
      y(1) = a(2).*b(3) - a(3).*b(2);
      y(2) = a(3).*b(1) - a(1).*b(3);
      y(3) = a(1).*b(2) - a(2).*b(1);
    end
    function y = dot(a,b)
      y = a(:)'*b(:);
    end
    function y = sum(a,dim)
      if length(size(a))>2
        error('sum works for 2d arrays only.')
      end
      if nargin==1
        if size(a,1) == 1
          dim = 2;
        else
          dim = 1;
        end
      end
      if dim==1
        y = a(1,:);
        for i=1:size(a,2)
          for j=2:size(a,1)
            y(i).c = y(i).c + a(j,i).c;
          end
        end
      else
        y = sum(a.',1).';
      end
    end
    function y = cumsum(a,dim)
      if length(size(a))>2
        error('cumsum works for 2d arrays only.')
      end
      if nargin==1
        if size(a,1) == 1
          dim = 2;
        else
          dim = 1;
        end
      end
      if dim==1
        y = a;
        for i=1:size(a,2)
          for j=2:size(a,1)
            y(j,i).c = y(j-1,i).c + a(j,i).c;
          end
        end
      else
        y = cumsum(a.',1).';
      end
    end
    function y = prod(a,dim)
      if length(size(a))>2
        error('prod works for 2d arrays only.')
      end
      if nargin==1
        if size(a,1) == 1
          dim = 2;
        else
          dim = 1;
        end
      end
      if dim==1
        y = a(1,:);
        for i=1:size(a,2)
          for j=2:size(a,1)
            y(i) = y(i).*a(j,i);
          end
        end
      else
        y = prod(a.',1).';
      end
    end
    function y = cumprod(a,dim)
      if length(size(a))>2
        error('cumprod works for 2d arrays only.')
      end
      if nargin==1
        if size(a,1) == 1
          dim = 2;
        else
          dim = 1;
        end
      end
      if dim==1
        y = a;
        for i=1:size(a,2)
          for j=2:size(a,1)
            y(j,i) = y(j-1,i).*a(j,i);
          end
        end
      else
        y = cumprod(a.',1).';
      end
    end
    function y = diff(a,n,dim)
      if length(size(a))>2
        error('diff works for 2d arrays only.')
      end
      if nargin<2
        n = 1;
      end
      if nargin<3
        if size(a,1) == 1
          dim = 2;
        else
          dim = 1;
        end
      end
      if dim == 1
        y = diff(eye(size(a,1)),n,dim)*a;
      else
        y = a*diff(eye(size(a,2)),n,dim);
      end
    end
    function y = diag(a,k)
      if nargin==1
        k = 0;
      end
      if isvector(a)
        I = diag(ones(size(a)),k);
        y = repmat(0*a(1),size(I));
        y(I>0) = a;
      else
        I = zeros(size(a));
        I(:) = 1:numel(I);
        y = a(diag(I,k));
      end
    end
    function a = tril(a,k)
      if nargin==1
        k = 0;
      end
      I = ones(size(a));
      a(~tril(I,k)) = 0*a(1);
    end
    function a = triu(a,k)
      if nargin==1
        k = 0;
      end
      I = ones(size(a));
      a(~triu(I,k)) = 0*a(1);
    end
    function d = det(A)
      [nA,mA] = size(A);
      if nA~=mA
        error('matrix must be square')
      elseif nA==1
        d = A;
        return
      elseif nA==2
        d = A(1).*A(4) - A(2).*A(3);
        return
      elseif nA==3
        d = (A([5 6 4]).*A([9 7 8])-A([6 4 5]).*A([8 9 7]))*A(:,1);
        return
      end
      d = 0.*A(1) + 1;
      while size(A,1)>3
        v = abs(aeval(A(1)));
        j = 1;
        if v(1)<1e-5
          A0 = aeval(aget(A,1));
          [~,j] = max(abs(A0(:)));
        end
        [i,m] = ind2sub(size(A),j);
        a = A(j);
        d = d.*(-1)^(i+m).*a;
        B = A(:,m)./a;
        C = A(i,:);
        B(i) = [];
        C(m) = [];
        A(i,:) = [];
        A(:,m) = [];
        A = A - B*C;
      end
      d = d.*det(A);
    end
    function [Q,A,E] = qr(A,mode)
      if nargout==3
        if nargin==1 || isequal(mode,'matrix')
          [~,~,E] = qr(aeval(aget(A,1)));
          A = A*E;
        elseif isequal(mode,'vector') || ~mode
          [~,~,E] = qr(aeval(aget(A,1)),'vector');
          A = A(:,E);
        else
          error('option not supported')
        end
      end
      if ~isreal(A)
        error('complex data not supported')
      end
      [nA,mA] = size(A);
      if nargin==1 || nA<=mA
        mode = 1;
      end
      if nargout>1
        Q = 0.*A(1) + eye(nA);
      end
      for i = 1:min(nA-1,mA)
        j1 = i:nA;
        j2 = i+1:mA;
        v = A(i:nA,i);
        z = -sign(v(1)+realmin)*norm(v);
        v = v./z;
        u1= sqrt(1-v(1));
        u = [u1;v(2:end)./(-conj(u1))];
        if nargout>1
          Q(:,j1) = Q(:,j1) - (Q(:,j1)*u)*u';
        else
          A(i+1:nA,i) = u(2:end)./u1;
        end
        A(i,i) = z;
        A(j1,j2) = A(j1,j2) - u*(u'*A(j1,j2));
      end
      if nargout>1
        A = triu(A);
      else
        Q = A;
      end
      if ~mode
        A = A(1:mA,1:mA);
        Q = Q(:,1:mA);
      end
    end
    function [L,R,P] = lu(A,mode)
      if diff(size(A))
        error('matrix must be square')
      end
      nA = size(A,1);
      A0 = aeval(A,0);
      [~,~,p] = lu(A0(:,:,1),'vector');
      A = A(p,:);
      L = 0*A(1) + eye(nA);
      R = repmat(0*A(1),nA,nA);
      for i = 1:nA
        R(i,i:nA) = A(i,i:nA) - L(i,1:i-1)*R(1:i-1,i:nA);
        L(i+1:nA,i) = (A(i+1:nA,i) - L(i+1:nA,1:i-1)*R(1:i-1,i))./R(i,i);
      end
      if nargout==2
        [~,q] = sort(p);
        L = L(q,:);
      elseif nargout==3
        if nargin==1 || ~isequal(mode,'vector')
          P = eye(nA);
          P = P(p,:);
        else
          P = p;
        end
      else
        error('nargout must be 2 or 3')
      end
    end
    function R = chol(A)
      nA = size(A,1);
      R = repmat(0*A(1),nA,nA);
      for i = 1:nA
        r = R(1:i-1,i);
        R(i,i) = sqrt(A(i,i) - r'*r);
        if ~isreal(R(i,i)) || R(i,i)==0
          error('matrix must be positive definite')
        end
        R(i,i+1:nA) = (A(i,i+1:nA) - r'*R(1:i-1,i+1:nA))./R(i,i);
      end
    end
    function y = inv(a)
      m = size(a,1);
      if size(a,2)~=m
        error('Matrix must be square.')
      end
      if m == 1
        y = 1./a;
      elseif m==2
        y = [a(4) -a(3);-a(2) a(1)]./(a(1).*a(4) - a(2).*a(3));
      else
        y = a\eye(m);
      end
    end
    function y = pinv(a)
      y = (a'*a)\a';
    end
    function [u,s,v] = svd(a)
      if diff(size(a))<0
        if nargout<2
          u = sqrt(real(eig(a'*a)));
        else
          [v,s,u] = svd(a');
          s = s';
        end
      else
        [u,s] = eig(a*a');
        s = sqrt(real(diag(s)))';
        if nargout<2
          u = s';
          return
        end
        v = (a'*u)./s;
        if diff(size(a))
          v = [v null(v')];
        end
        s = [diag(s) zeros(size(a,1),diff(size(a)))];
      end
    end
    function varargout = eig(a,b)
      if nargin==2
        varargout = cell(1,nargout);
        [varargout{:}] = eig(b\a);
      else
        na = size(a,1);
        if na==1
          if nargout==1
            varargout = {a};
          else
            varargout = {1,a};
          end
        elseif na==2
          p = trace(a)/2;
          d = sqrt(p.*p-det(a));
          d = [p+d;p-d];
          if nargout<2
            varargout = {d};
          else
            for i = 2:-1:1
              aa = a - d(i)*eye(2);
              [~,~,e] = qr(aeval(aget(aa,1))',0);
              v(:,i) = [-aa(e(1),2); aa(e(1),1)];
            end
            v = v./v(2,:);
            v = v./[norm(v(:,1)) norm(v(:,2))];
            varargout = {v,diag(d)};
          end
        elseif na==3
          q = a([2 3 6])*a([4 7 8]).' - a([1 5 9])*a([5 9 1]).';
          d = roots([-1 trace(a) q det(a)]);
          if ishermitian(a)
            d = real(d);
          end
          if nargout==1
            varargout = {d};
          else
            for i = 3:-1:1
              aa = a - d(i)*eye(3);
              [~,~,e] = qr(aeval(aget(aa,1))',0);
              v(i,:) = cross(aa(e(1),:),aa(e(2),:));
            end
            v = (v./v(:,3)).';
            v = v./sqrt(sum(conj(v).*v,1));
            varargout = {v,diag(d)};
          end
        else
          e0 = eye(na);
          a0 = aeval(a,0);
          [v0,d0] = eig(a0);
          a = v0\(a*v0);
          da= a - aeval(a,0);
          d = repmat(0.*a(1),na,1);
          v = repmat(0.*a(1),na,na);
          for i = 1:na
            l = d0(i,i);
            w = e0(:,i);
            q = 1./(diag(d0)-l);
            q(i) = 0;
            for j = 1:a(1).k
              l = d0(i,i) + (da(i,:)*w)./w(i);
              w = w + q.*(l*w-a*w);
            end
            d(i) = l;
            w = v0*w;
            v(:,i) = w./norm(w);
          end
          if nargout==1
            varargout = {d};
          else
            varargout = {v,diag(d)};
          end
        end
      end
    end
    function y = orth(a)
      [y,~] = qr(a,0);
    end
    function y = null(a)
      [na,ma] = size(a);
      if na>=ma
        y = a([]);
      else
        [~,R] = qr(a,0);
        y = orth([-R(:,1:na)\R(:,na+1:ma);0.*a(1)+eye(ma-na)]);
      end
    end
    function y = trace(a)
      y = 0.*a(1);
      for i = 1:size(a,1)
        y = y + a(i,i);
      end
    end
    function y = norm(a,p)
      if isempty(a)
        y = a;
        return
      end
      if numel(a)==1
        y = abs(a);
        return
      end
      if nargin==1
        p = 2;
      end
      if ~isvector(a)
        if isinf(p)
          y = max(sum(abs(a),2));    % maximal row sum
        elseif p==1
          y = max(sum(abs(a),1));    % mximal column sum
        elseif p==2
          if size(a,1)<size(a,2)
            y = sqrt(max(eig(a*a')));
          else
            y = sqrt(max(eig(a'*a)));
          end
        elseif isequal(p,'fro')
          y = norm(a(:),2);
        else
          error('The only matrix norms available are 1, 2, inf, and ''fro''.')
        end
      else
        if isinf(p)
          if p>0
            y = max(abs(a));
          else
            y = min(abs(a));
          end
        elseif p==2
          y = sum(conj(a).*a).^(1/2);
        else
          y = sum(abs(a).^p).^(1/p);
        end
      end
    end
    function y = expm(a)
      [v,d] = eig(a);
      y = v.*(exp(diag(d)).');
      y = y/v;
    end
    function y = logm(a)
      [v,d] = eig(a);
      y = v.*(log(diag(d)).');
      y = y/v;
    end
    function y = sqrtm(a)
      [v,d] = eig(a);
      y = v.*(sqrt(diag(d)).');
      y = y/v;
    end
    %% BINARY OPERATOS
    function y = plus(a,b)
      if isa(b,'audi') && ~isempty(b) && isequaln(nan,b.s) % special case stand-alone
        y = b;
        y.c = b.c(ones(numel(a),1),:);
        y.c(:,1) = a(:);
        y.s = size(a);
        return
      end
      if ~isa(a,'audi')
        y = plus(b,a);
        return
      end
      % now, a is audi
      if isempty(a) || isempty(b)
        y = a([]);
        return
      end
      [a,b] = singexp(a,b);
      if ~isequal(size(a),size(b))
        error('dimension mismatch')
      end
      if ~isa(b,'audi')
        y = a;
        for i = 1:numel(y)
          y(i).c(:,1) = a(i).c(:,1) + b(i);
        end
      else    % both a and b are audis
        if a(1).k <= b(1).k
          j = 1:size(a(1).c,2);
          y = a;
        else
          j = 1:size(b(1).c,2);
          y = b;
        end
        % now, a and b have equal order
        for i = 1:numel(a)
          y(i).c = a(i).c(:,j) + b(i).c(:,j);
        end
      end
    end
    function y = minus(a,b)
      y = a + (-b);
    end
    function [y,yc] = times(a,b)
      if isempty(a) || isempty(b)
        y = a([]) + b([]);
        yc = [];
        return
      end
      [a,b] = singexp(a,b);
      yc = 0;           % data matrix if y is audi array of scalars
      if ~isa(a,'audi') || (isa(b,'audi') && numel(b)>numel(a))
        [y,yc] = times(b,a);
        return
      end
      if any(size(a)~=size(b))
        error('cannot multiply arrays of different size')
      end
      % now, a is audi, and a and b have equal size
      if ~isa(b,'audi')              % audi times numeric
        y = a;
        for i = 1:numel(a)
          y(i).c = y(i).c*b(i);
        end
      else                           % audi times audi
        if a(1).k < b(1).k
          j = 1:size(a(1).c,2);
          for i = 1:numel(a)
            b(i).c = b(i).c(:,j);
            b(i).k = a(1).k;
          end
          y = a;
        elseif b(1).k < a(1).k
          j = 1:size(b(1).c,2);
          for i = 1:numel(a)
            a(i).c = a(i).c(:,j);
            a(i).k = b(1).k;
          end
          y = b;
        else
          y = a;
        end
        % now, a and b have equal order
        [~,F] = audi.mat(a(1).k,a(1).n);
        I = sumidx(a(1));
        iF= 1./F;
        if numel(a(1).c(:,1))==1                % asize==[1 1]
          A = zeros(numel(a),size(a(1).c,2));
          B = A;
          for i = 1:numel(a)
            A(i,:) = a(i).c;
            B(i,:) = b(i).c;
          end
          A = iF.*A;
          B = iF.*B;
          Y = A;
          for j = 1:size(y(1).c,2)
            Y(:,j) = sum(A(:,I{j}(:,1)).*B(:,I{j}(:,2)),2);
          end
          yc = F.*Y;
          for i = 1:numel(a)
            y(i).c = yc(i,:);
          end
        else                                   % asize~=[1 1]
          for i = 1:numel(y)
            ac = iF.*a(i).c;
            bc = iF.*b(i).c;
            yc = ac;
            for j = 1:size(y(1).c,2)
              yc(:,j) = sum(ac(:,I{j}(:,1)).*bc(:,I{j}(:,2)),2);
            end
            y(i).c = F.*yc;
          end
        end
      end
    end
    function y = rdivide(a,b)
      y = a.*(b.^(-1));
      if numel(a)*numel(b)==1 && isa(a,'audi') && isa(b,'audi') && adim(a)*adim(b)==1
        i = a.c(:,1)==0 & b.c(:,1)==0;
        if any(i)
          aa = aget(a,i);
          bb = aget(b,i);
          aa.c = [aa.c(:,2:end)./(1:aa.k) nan(size(aa.c,1),1)];
          bb.c = [bb.c(:,2:end)./(1:bb.k) nan(size(bb.c,1),1)];
          y = aset(y,aa./bb,i);
        end
      end
    end
    function y = ldivide(a,b)
      y = rdivide(b,a);
    end
    function y = mtimes(a,b)
      if isempty(a) || isempty(b)
        if ~size(a,2) && ~size(b,1)
          y = zeros(size(a,1),size(b,2));
        else
          y = a([]) + b([]);
        end
      elseif numel(a)==1 || numel(b)==1
        % multiplication by scalar
        y = a.*b;
      elseif isnumeric(b)
        % only first argument can be numeric
        y = (mtimes(b.',a.')).';
      elseif size(a,2)==size(b,1)
        % matrix multiplication
        y = repmat(0*b(1),size(a,1),size(b,2)); % initialize result
        if size(a,2)==1 && size(b,1)==1         % column times row
          if size(a,1)<size(b,2)
            for i = 1:size(a,1)
              y(i,:) = a(i).*b;
            end
          else
            for i = 1:size(b,2)
              y(:,i) = a.*b(i);
            end
          end
        elseif isnumeric(a) && size(b(1).c,1)==1
          % numerical array times audi array of scalars
          B = zeros(size(b,1),size(b(1).c,2));
          for j=1:size(b,2)
            for i=1:size(b,1)
              B(i,:) = b(i,j).c;
            end
            BB = a*B;
            for i=1:size(a,1)
              y(i,j).c = BB(i,:);
            end
          end
        else
          for i = 1:size(a,1)
            A = a(i,:).';
            for j = 1:size(b,2)
              if size(b(1).c,1)==1
                [yy,yc] = times(A,b(:,j));
                y(i,j) = yy(1);
                y(i,j).c = sum(yc,1);
              else
                y(i,j) = sum(A.*b(:,j));
              end
            end
          end
        end
      else
        error('Inner matrix dimensions must agree.')
      end
    end
    function y = mldivide(a,b)
      if isnumeric(a)
        y = (a\eye(size(a,1)))*b;
        return
      elseif numel(a)==1
        y = ldivide(a,b);
        return
      elseif size(a,1) ~= size(b,1)
        error('dimension mismatch')
      end
      % now, a and b are audis of compatible size
      if size(a,1)==size(a,2)     % square
        if isdiag(a)
          y = diag(a).\b;
        elseif istril(a)
          y = repmat(0.*a(1),size(a,1),size(b,2));
          for i = 1:size(a,1)
            y(i,:) = (b(i,:) - a(i,1:i-1)*y(1:i-1,:))./a(i,i);
          end
        elseif istriu(a)
          y = a(end:-1:1,end:-1:1)\b(end:-1:1,:);
          y = y(end:-1:1,:);
        else
          [L,U,p] = lu(a,'vector');
          y = U\(L\b(p,:));
        end
      elseif size(a,1)>=size(a,2) % overdetermined
        [Q,R] = qr(a);
        b = Q'*b;
        y = R(1:size(a,2),:)\b(1:size(a,2),:);
      else                        % underdetermined
        a0= aeval(a,0);
        f = any(a0(:,:,1)\eye(size(b,1)),2);
        y = repmat(0.*a(1),size(a,2),size(b,2));
        y(f,:) = a(:,f)\b;
      end
    end
    function y = mrdivide(a,b)
      y = mldivide(b.',a.').';
    end
    function y = power(a,b)
      if isnumeric(b)
        if b==0
          y = 0.*a + 1;
        elseif b==1
          y = a;
        elseif b==2
          y = a.*a;
        elseif b==3
          y = a.*a.*a;
        elseif b==4
          y = a.*a;
          y = y.*y;
        elseif b==-1
          y = repmat(0.*a(1),size(a));
          for i=1:numel(a)
            h = 1./a(i).c(:,1);
            x =-a(i);
            x.c(:,1) = 0;
            h = h.^(a(i).k+1);
            y(i).c(:,1) = h;
            for j = 1:a(i).k
              h = h.*a(i).c(:,1);
              y(i) = y(i).*x;
              y(i).c(:,1) = y(i).c(:,1) + h;
            end
          end
        elseif b==-2
          y = power(a,-1);
          y = y.*y;
        elseif b==-3
          y = power(a,-1);
          y = y.*y.*y;
        elseif b==-4
          y = power(a,-1);
          y = y.*y;
          y = y.*y;
        else
          df = @(x,m) m.*power(x,m-1);
          y = adval(df,a,b,1);
        end
      else
        y = exp(b.*log(a));
      end
      %       [a,b] = singexp(a,b);
      %       for i = 1:numel(y)
      %         y(i).c(:,1) = aeval(a(i)).^aeval(b(i));
      %       end
    end
    function y = mpower(a,b)
      if isnumeric(b) && b==round(b)
        if b<0
          y = mpower(inv(a),-b);
        elseif b==0
          y = 0*a + eye(size(a));
        elseif b==1
          y = a;
        elseif b==2
          y = a*a;
        elseif b<8
          y = a*mpower(a,b-1);
        else
          bb = floor(b/2);
          y = mpower(a,b-2*bb)*mpower(mpower(a,bb),2);
        end
      else
        error('Exponent must be positive integer.')
      end
    end
    function [y,i] = min(a,b,dim)
      % MIN Smallest component.
      %    For univariate real a and b, derivatives of min(a,b)
      %    are taken from the right.
      %
      %    Example: If x = ainit(0,2) then
      %                aeval(min(0,x),1) = 0
      %                aeval(min(0,-x^2),2) = -2
      %    See also: audi.max
      if nargin==1
        if isvector(a)
          [~,i] = min(aeval(a,0));
          y = a(i);
        else
          [y,i] = min(a,[],1);
        end
      elseif nargin==2
        i = [];
        if numel(a)==1 && numel(b)==1
          [a,b] = adapt(a,b);
          a0 = aeval(a,0);
          b0 = aeval(b,0);
          y = a;
          if isreal(a) && isreal(b)
            if adim(a)==1
              j = a0>b0;
              e = a0==b0;
              for d = 1:aord(a)
                j = j | (e & (aeval(a,d)>aeval(b,d)));
                e = e & (aeval(a,d)==aeval(b,d));
              end
            else
              j = a0>b0;
            end
          else
            j = abs(a0)>abs(b0)
          end
          y = aset(y,aget(b,j),j);
        else
          [a,b] = singexp(a,b);
          y = 0.*a;
          for j = 1:numel(a)
            y(j) = min(a(j),b(j));
          end
        end
      else
        if length(size(a))>3
          error('2d arrays only')
        end
        [~,i] = min(aeval(a,0),[],dim);
        y = repmat(0.*a(1),size(i));
        if dim==2
          a = a';
        end
        for j = 1:length(i)
          y(j) = a(i(j),j);
        end
      end
    end
    function [y,i] = max(a,b,c)
      % MAX Largest component.
      %    For univariate real a and b, derivatives of min(a,b)
      %    are taken from the right.
      %
      %    Example: If x = ainit(0,2) then
      %                aeval(max(0,x^2),2) = 2
      %                aeval(max(0,-x),1) = 0
      %    See also: audi.min
      if nargin==1
        [y,i] = min(-a);
      elseif nargin==2
        [y,i] = min(-a,-b);
      else
        [y,i] = min(-a,b,-c);
      end
      y = -y;
    end
    %% COMPLEX NUMBERS
    function y = real(a)
      y = a;
      for i = 1:numel(a)
        y(i).c = real(a(i).c);
      end
    end
    function y = imag(a)
      y = a;
      for i = 1:numel(a)
        y(i).c = imag(a(i).c);
      end
    end
    function y = conj(a)
      y = a;
      for i = 1:numel(a)
        y(i).c = conj(a(i).c);
      end
    end
    function y = angle(a)
      y = imag(log(a));
    end
    %% CHECKING PROPERTIES
    function y = isreal(a)
      y = true;
      for i = 1:numel(a)
        y = y && isreal(a(i).c);
      end
    end
    function y = isfinite(a)
      y = true(size(a));
      for i = 1:numel(a)
        y(i) = all(isfinite(a(i).c(:)));
      end
    end
    function y = istril(a)
      f = find(triu(true(size(a)),1));
      y = true;
      for i = 1:numel(f)
        if any(a(f(i)).c(:))
          y = false;
          return
        end
      end
    end
    function y = istriu(a)
      y = istril(a.');
    end
    function y = isdiag(a)
      y = istril(a) && istril(a.');
    end
    function y = issymmetric(a)
      y = isequal(a - a.',0.*a);
    end
    function y = ishermitian(a)
      y = isequal(a - a',0.*a);
    end
    %% QUANTIFIERS
    function y = all(a,varargin)
      y = all(a~=0,varargin{:});
    end
    function y = any(a,varargin)
      y = any(a~=0,varargin{:});
    end
    %% FOURIER
    function x = fft(x,varargin)
      for i = 1:size(x.c,2)
        x.c(:,i) = fft(x.c(:,i),varargin{:});
      end
    end
    function x = ifft(x,varargin)
      for i = 1:size(x.c,2)
        x.c(:,i) = ifft(x.c(:,i),varargin{:});
      end
    end
    %% UNARY OPERATORS
    function y = uminus(a)
      y = a;
      for i = 1:numel(a)
        y(i).c = -y(i).c;
      end
    end
    function y = uplus(a)
      y = a;
    end
    function y = ctranspose(a)
      y = builtin('transpose',conj(a));
    end
    %% RELATIONAL OPERATORS
    function y = le(a,b)
      y = a<b | a==b;
    end
    function y = ge(a,b)
      y = a>b | a==b
    end
    function y = lt(a,b)
      y = lt(aeval(a,0),aeval(b,0));
    end
    function y = gt(a,b)
      y = gt(aeval(a,0),aeval(b,0));
    end
    function y = eq(a,b)
      if isnumeric(b)
        b = b + 0*a(1);
      end
      if isnumeric(a)
        a = a + 0*b(1);
      end
      if numel(a)==1
        a = repmat(a,size(b));
      end
      if numel(b)==1
        b = repmat(b,size(a));
      end
      y = true(size(a));
      for i=1:numel(a)
        y(i) = isequal(a(i),b(i));
      end
    end
    function y = ne(a,b)
      y = ~(a==b)
    end
    function y = and(a,b)
      y = and(aeval(a,0),aeval(b,0));
    end
    function y = or(a,b)
      y = or(aeval(a,0),aeval(b,0));
    end
    function y = xor(a,b)
      y = xor(aeval(a,0),aeval(b,0));
    end
    function y = not(a)
      y = not(aeval(a,0));
    end
    %% CATENATION AND EVALUATION
    function y = horzcat(varargin)
      v = audi.const2audi(varargin);
      y = builtin('horzcat',v{:});
    end
    function y = vertcat(varargin)
      v = audi.const2audi(varargin);
      y = builtin('vertcat',v{:});
    end
    function y = cat(dim,varargin)
      v = audi.const2audi(varargin);
      y = builtin('cat',dim,v{:});
    end
    function y = subsref(a,s)
      if ~strcmp(s(1).type,'{}')
        y = builtin('subsref',a,s);
      else
        if isempty(s(1).subs{end})
          y = aeval(a,[s(1).subs{1:end-1}],[]);
        else
          y = aeval(a,[s(1).subs{:}]);
        end
        if numel(s)>1
          y = builtin('subsref',y,s(2:end));
        end
      end
    end
    %% AUXILIARY ROUTINES
    function [aa,bb] = singexp(a,b)
      if isequal(size(a),size(b))
        aa = a;
        bb = b;
      else
        aa = repmat(a,(size(a)~=1)+(size(a)==1).*size(b));
        bb = repmat(b,(size(b)~=1)+(size(b)==1).*size(a));
      end
    end
    function c = bsxfun(f,a,b)
      c = f(repmat(a,(size(a)~=1)+(size(a)==1).*size(b)),...
        repmat(b,(size(b)~=1)+(size(b)==1).*size(a)));
    end
    function varargout = intern_taylor(D,x0,y0)
      if nargin==1
        a = D;
        if prod(asize(a)) > 1
          error('Evaluation at single points only.')
        elseif size(a,2) > 1
          error('Evaluation for column vectors only.')
        elseif adim(a) > 1
          error('Evaluation for univariate functions only.')
        else
          y = zeros(size(a,1),aord(a)+1);
          for i = 1:size(a,1)
            b = a(i);
            for j = 0:b.k
              y(i,end-j) = b.c(:,j+1)/factorial(j);
            end
          end
          [i,j] = find(isnan(y));
          if ~isempty(i)
            y = y(:,max(j)+1:end);
          end
        end
        varargout{1} = y;
      else
        if ~iscell(D)
          D = {D};
        end
        N = numel(D);
        for i = 1:N
          z = D{i};
          for j = 1:N
            yold = ['y' num2str(j)];
            ynew = ['audi.odeaux(x,' num2str(j) ',D,y0)'];
            z = strrep(z,yold,ynew);
          end
          %           if N==1
          %             s = strrep(s,'y','audi.odeaux(x,1,D,y0)');
          %           end
          D{i} = str2func(['@(x,b,D,y0) ' z]);
        end
        varargout = cell(1,min(N,nargout));
        for i = 1:min(N,nargout)
          varargout{i} = intern_taylor(audi.odeaux(x0,i,D,y0));
        end
      end
    end
    function i = idx(a,s)
      if a.k <= 2
        if sum(s)==0
          i = 1;
        elseif sum(s)==1
          i = 1 + find(s);
        elseif max(s)==2
          i = 1 + a.n + find(s);
        else
          f = find(s);
          i = f(2) - f(1)*(f(1)+1)/2 + a.n*(f(1)+1) + 1;
        end
      else
        J = audi.mat(a.k,a.n);
        i = find(all(J==s,2));
      end
    end
    function s = sub(a,i)
      J = audi.mat(a.k,a.n);
      s = J(i,:);
    end
    function I = sumidx(a)
      global II
      kk = a.k;
      nn = a.n;
      if kk==0
        I = {[1 1]};
        return
      end
      try
        I = II{kk,nn};
        if isempty(I)
          error('not defined yet')
        end
      catch
        J = audi.mat(kk,nn);
        I = cell(1,size(J,1));
        I{1} = [1 1];
        for i = 2:size(J,1)
          p = J(i,:);
          f = [find(all(J(1:nchoosek(nn+sum(p)-1,nn),:)<=p,2));i];
          I{i} = [f f];
          for j = 1:numel(f)
            I{i}(j,2) = idx(a,p-J(f(j),:));
          end
        end
        II{kk,nn} = I;
      end
    end
    function [a,b] = adapt(a,b)
      %ADAPT Generate compatible audi scalars.
      %   [aa,bb] = ADAPT(a,b) returns two compatible audis.
      %   Input arguments of class double become audi constants.
      %   If a and b are audi scalars, aa and bb have equal order,
      %   aord(aa) = aord(bb) = min(aord(a),aord(b)).
      if isempty(a) || isempty(b)
        if isa(a,'audi')
          a = a([]);
          b = a;
        else
          b = b([]);
          a = b;
        end
      elseif isnumeric(a)
        a = 0*b + a;
      elseif isnumeric(b)
        b = 0*a + b;
      elseif a.n ~= b.n
        error('Audi variables must have equal dimension.')
      elseif a.k < b.k
        aa = a;
        for i = 1:size(a.c,2)
          aa.c(:,i) = b.c(:,idx(b,sub(aa,i)));
        end
        b = aa;
      elseif a.k > b.k
        [b,a] = adapt(b,a);
      end
    end
    function [a1,a2,p] = prepfct(~,d,k)
      if nargin==2
        k = sum(d);
      else
        d = double(d);
        d(d==false) = nan;
      end
      i = isfinite(d);
      j = isnan(d);
      nd= numel(d);
      p = cell(1,nd);
      [p{i}] = deal(0);
      [p{i}] = ainit(p{i},k);
      [p{j}] = deal(0);
      %       d = d(i);
      a1 = '';
      a2 = '';
      for i = 1:nd
        a1 = [a1 ',x' num2str(i)];
        a2 = [a2 ',x' num2str(i) '+p{' num2str(i) '}'];
        p{i}.s = nan;       % used to detect special case when comuting xi+p{i}
      end
      a1 = a1(2:end);
      a2 = a2(2:end);
    end
  end
  %% Static
  methods(Static)
    function a = const2audi(a)
      for i = 1:numel(a)
        if ~isempty(a{i}) && ~isnumeric(a{i})
          z = 0*a{i}(1);
          break
        end
      end
      for i = 1:numel(a)
        a{i} = a{i} + z(ones(size(a{i})));
      end
    end
    function [J,F] = mat(k,n)
      persistent JJ FF
      if k==0
        J = sparse(1,n);
        F = 1;
        return
      end
      try
        J = JJ{k,n};
        F = FF{k,n};
        if isempty(J) || isempty(F)
          error('not defined yet')
        end
      catch
        % J and F not computed yet
        if k==1
          J = [sparse(1,n);speye(n)];
        elseif k==2
          J = [sparse(1,n);speye(n);2*speye(n)];
          for i = 1:n-1
            J = [J;sparse(n-i,i-1) sparse(n-i,1)+1 speye(n-i)];
          end
        else
          v = cell(1,n);
          [v{:}] = ndgrid(0:k);
          V = zeros(numel(v{1}),n);
          for i = 1:n
            V(:,i) = v{i}(:);
          end
          J = audi.mat(2,n);
          for kk = 3:k
            J = [J;V(sum(V,2)==kk,:)];
          end
        end
        F = full(prod(factorial(J),2))';
        JJ{k,n} = J;
        FF{k,n} = F;
      end
    end
    function f = odeaux(t,b,D,y0)
      if isnumeric(t)
        f = y0(b);
      else
        f  = adval(D{b},t,b,D,y0,1);
      end
    end
  end
end