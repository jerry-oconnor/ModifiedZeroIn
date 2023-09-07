% -------------------------------------------------------------------------
% Modified ZeroIn
%   modification of brent's method, utilizes a combination of 
%   bisection and IQI 
%
%   Based on the method defined in the following paper: 
%       Wilkins, Gautam & Gu, Ming. (2013). 
%       A modified Brent’s method for finding zeros of functions. 
%       Numerische Mathematik. 123. 10.1007/s00211-012-0480-x. 
%
%       Link: 
%           https://dl.acm.org/doi/abs/10.1007/s00211-012-0480-x
% -------------------------------------------------------------------------
% Inputs: 
%   f : f(x) with x ∈ ℝ
%   int : object that denotes an interval [int.a, int.b]
%   params: object with the following fields: 
%               params.func_tol: function tolerance 
%               params.root_tol: root tolerance 
%               params.maxit: maximum number of iterations 
% -------------------------------------------------------------------------
% Outputs: 
%   root : computed root or zero 
%   info.it : number of iterations 
%   info.flag : success flag, 0 if successful, 1 if unsuccessful 
% -------------------------------------------------------------------------

% Main function 
function [root,info] = modifiedzeroin(func, Int, params)
    % initial vals 
    a = Int.a; b = Int.b;
    info.it = 0; 
    
    xHist = [];
    m = containers.Map('KeyType','double', 'ValueType','double');

    % begin iteration 
    x0 = a; x1 = b; x2 = ((a + b) / 2); x3 = x2; 
    while (abs(x1 - x0) > params.root_tol)                ...
            && abs(mapCheck(m,func,x3)) > params.func_tol ...
            && info.it < params.maxit

        % calc x3 via iqi 
        x3 = IQI(func,x0,x1,x2,m);
        
        % check inclusion of x3 
        if x3 < x0 || x3 > x1 % x3 ∉ [a,b]
            [x1, x2, x3] = bisection(func,x0,x1,params.root_tol,7,m);
        end
        xHist = [xHist, x3]; % add new x3 to list 

        % update x values and interval 
        x0 = x1; x1 = x2; x2 = x3;         
        a_n = min([x0,x1,x2]); b_n = max([x0,x1,x2]); 
        x0 = a_n; x1 = b_n; x2 = (x0 + x1) / 2; 
                
        % bisect for f3 if not decreasing by at least factor of 2 
        decDist = 4; 
        if rem(info.it - 1, decDist) == 0 && info.it > 1
            x3old = xHist(info.it - decDist);

            f3 = mapCheck(m,func,x3); 
            f3old = mapCheck(m,func,x3old);
            
            % decrease by a factor of 2
            if abs(f3old / f3) < 2
                [x1, x2, x3] = bisection(func,x0,x1,params.root_tol, 3, m);
            end
        end 
        info.it = info.it + 1;
    end
    
    % determine failure or success 
    if info.it == params.maxit
        info.flag = 1;
    else
        info.flag = 0;
    end

    % return computed root 
    root = x3; 
    return; 
end 


% *************************************************************************
%  _    _  _    _  _   ______                    _    _                     
% | |  | || |  (_)| | |  ____|                  | |  (_)                    
% | |  | || |_  _ | | | |__  _   _  _ __    ___ | |_  _   ___   _ __   ___  
% | |  | || __|| || | |  __|| | | || '_ \  / __|| __|| | / _ \ | '_ \ / __| 
% | |__| || |_ | || | | |   | |_| || | | || (__ | |_ | || (_) || | | |\__ \ 
%  \____/  \__||_||_| |_|    \__,_||_| |_| \___| \__||_| \___/ |_| |_||___/ 
% *************************************************************************
                                                                             
% mapcheck 
%   - checks if x : func(x) is in map m, returns m(x) if yes, f(x) if not
%   - adds x : func(x) to m if not found 
function f = mapCheck(m, func, x)
    if isKey(m,x)
        f = m(x);
    else 
        f = func(x);
        m(x) = f;
    end
end 

% IQI 
%   - calculates inverse quadratic interpolation for x0-2
function x3 = IQI(f,x0,x1,x2,m)
    
    % function values 
    f0 = mapCheck(m, f, x0); 
    f1 = mapCheck(m, f, x1); 
    f2 = mapCheck(m, f, x2);
    
    % calculate and return x3 
    x3 = (((f1*f2) / ((f0 - f1)*(f0 - f2)))*x0) ...
         + (((f0*f2) / ((f1 - f0)*(f1 - f2)))*x1) ...
         + (((f0*f1) / ((f2 - f0)*(f2 - f1)))*x2);  
end

% Bisection 
%   - performs bisection for given interval, tol, and number of iterations 
function [x1, x2, x3] = bisection(f,a,b,tol,maxIt,m)
    err = 1;
    i = 0;
    while (i <= maxIt && err >= tol) 
        i = i + 1;
        p = (b + a)/2; 
        fa = mapCheck(m, f, a);
        fp = mapCheck(m, f, p);  
        error = abs(fp); 
        %  error = abs(x - xold); 
        if (error < tol || i >= maxIt) 
            %x3 = p; 
            return; 
        else 
            if (fa*fp < 0) 
                % root is between a and x 
                b = p; 
            else 
                %  root is between x and b 
                a = p; 
            end 
        end 
        x1 = a; x2 = b; x3 = (a + b) / 2;
    end
end 




