function K = presekDvehPloskev (f1, C1, f2, C2, gradf1, gradf2, X0, korak, funkcija)
  % Funkcija presehDvehPloskev (f1, C1, f2, C2, gradf1, gradf2, X0, korak, funkcija) za vhodne podatke:
  % f1       ... Funkcija prve ploskve
  % C1       ... Vrednost prve funkcije
  % f2       ... Funkcija druge ploskve
  % C2       ... Vrednost druge funkcije
  % gradf1   ... Gradient prve ploskve - podan kot stolpicni vektor
  % gradf1   ... Gradient prve ploskve - podan kot stolpicni vektor
  % X0       ... Zacetni priblizek
  % korak    ... Korak pri Eulerjevi oz RK4 metodi
  % funkcija ... Ce 1: Euler
  %              Ce 0: RK4
  % izracuna tocke na krivulji preseka in jih vrne v vektorju K
  
  % F(x) = (graf f1(x)) x (grad f2(x)))  /  (|| (graf f1(x)) x (grad f2(x))) ||)
  crossGrad = @(x) cross(gradf1(x),gradf2(x));
  F = @(x)  crossGrad(x) / norm(crossGrad(x));  % (|| (graf f1(x)) x (grad f2(x))) ||) == norm(crossGrad(x)) ??
  
  y = 0;% <-- priblizek
  v = F(y);
  % Jacobijeva matrika
  %       |  grad f1  |
  %   JG =|  grad f2  |
  %       | grad (v*x)| <- v' (v = F(y))
  JG = [gradf1; gradf2; v'];
  
  % IZRACUN TOCK KRIVULJE PRESEKA
  K = NaN;
  
  % Izracun po Eulerjevi metodi
  if(funkcija == 1)
    % euler(f ( dY/dt = f(t, Y)), interval [t0, tk], zacetni pogoj (Y(t0) = Y0), st. korakov)
    %[t, Y] = euler(f, tspan, Y0, n);
  
  % Izracun po RK4
  elseif(funkcija == 0)
  
  
  endif
  
  % Narisemo obe funkciji ter presecisce
  [y,x,z] = ndgrid(linspace(-5,5,64));
  cla reset;
  f = f1(x,y,z);
  g = f2(x,y,z);
  isosurface(x,y,z,f,C1);
  hold on;
  isosurface(x,y,z,g,C2);
  axis equal;
endfunction

function [t, Y] = euler(f, tspan, Y0, n)
  % [t, Y] = euler(f, tspan, Y0, n) poisce priblizno resitev DE 
  % dY/dt = f(t, Y) z zacetnim pogojem Y(t0) = Y0 
  % z Eulerjevo metodo z n koraki na intervalu tspan = [t0, tk].

  h = (tspan(2) - tspan(1))/(n-1);
  t = linspace(tspan(1), tspan(2), n);
  Y = Y0;
  for k = 1:(n-1)
    Y(:,k+1) = Y(:,k) + feval(f, t(k), Y(:,k))*h;
  end
endfunction

function [t, Y] = rk4(f, tspan, Y0, n)
  % [t, Y] = rk4(f, tspan, Y0, n) poisce priblizno resitev DE 
  % dY/dt = f(t, Y) z zacetnim pogojem Y(t0) = Y0 
  % z Runge-Kutta metodo 4. reda z n koraki na intervalu 
  % tspan = [t0, tk].

  h = (tspan(2) - tspan(1))/(n-1);
  t = linspace(tspan(1), tspan(2), n);
  Y = Y0;
  for k = 1:(n-1)
    k1 = h*feval(f, t(k), Y(:,k));
    k2 = h*feval(f, t(k) + h/2, Y(:,k) + k1/2);
    k3 = h*feval(f, t(k) + h/2, Y(:,k) + k2/2);
    k4 = h*feval(f, t(k) + h, Y(:,k) + k3);
    Y(:,k+1) = Y(:,k) + (k1 + 2*k2 + 2*k3 + k4)/6;
  endfor
endfunction

function [x, k] = newton(F, JF, x0, tol, maxit)
  %[x, k] = newton(F, JF, x0, tol, maxit) poisce priblizek 
  %x za resitev enacbe F(x) = 0 z Newtonovo metodo.
  %(k je stevilo korakov, ki jih metoda porabi.)
  %F... funkcija
  %JF... Jacobijeva matrika F
  %x0... zacetni priblizek
  %tol... zahtevana natancnost
  %maxit... najvecje stevilo iteracij

  for k = 1:maxit
    %Izvedemo en korak Newtonove metode...
    x = x0 - feval(JF, x0)\feval(F, x0);
    %... in testiramo, ce je metoda ze 'konvergirala'.
    if(norm(x - x0) < tol)
      break;
    end
    x0 = x;
  end

  %Izpisemo opozorilo v primeru, da zadnji priblizek ni znotraj tolerancnega obmocja.
  if(k == maxit)
    disp("Warning: The method did not converge after maxit iterations.")
  end
endfunction