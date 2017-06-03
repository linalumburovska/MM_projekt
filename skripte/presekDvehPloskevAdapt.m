function K = presekDvehPloskevAdapt(f1, C1, f2, C2, gradf1, gradf2, X0, korak, funkcija)
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
  
  % ZACASNO: vhodni podatki
  X0=[1;0.5;1];
     f1 = @(a) e.^(-a(1).^2+1) + a(2).^2+a(3).^ 2;
   f2 = @(a) a(1).^2 + a(2).^2 + a(3).^2;

   gradf1 =@(a) [-2*a(1)*e.^(-a(1).^2+1);2*a(2);2*a(3)];
   gradf2 =@(a) [2*a(1);2*a(2);2*a(3)];



  % F(x) = (graf f1(x)) x (grad f2(x)))  /  (|| (graf f1(x)) x (grad f2(x))) ||)
  crossGrad = @(x) cross(gradf1(x),gradf2(x));
  F = @(x)  crossGrad(x) / norm(crossGrad(x));

  % Jacobijeva matrika
  %       |  grad f1  |
  %   JG =|  grad f2  |
  %       | grad (v*x)| <- v' (v = F(y))
  JG = @(y)[gradf1(y)'; gradf2(y)'; F(y)'];
  
  G = @(X) [f1(X)-C1; f2(X)-C2; F(Y)'*X-F(Y)'*Y];
  
  % IZRACUN TOCK KRIVULJE PRESEKA
  K = NaN;
  
  % postavitev zacetnega priblizka na presek
  v = F(X0); 
  JGv = @(X)[gradf1(X)'; gradf2(X)'; v'];
  Gv = @(X) [f1(X)-C1; f2(X)-C2; v'*X-v'*X0];
  X0 = newton(Gv, JGv, X0, 10^-8, 20); 

  % Izracun po Eulerjevi metodi
  if(funkcija == 1)
    % euler(f ( dY/dt = f(t, Y)), interval [t0, tk], zacetni pogoj (Y(t0) = Y0), st. korakov)
    % y' = f (t, y) <-- enacba
    % y(t0) = y0    <-- zacetni pogoj
    % Y             <-- nabor funkcijskih vrednosti
    % t             <-- ob casih t
    f = @(t,Y) [F(Y)];
    [t, Y] = euler(f, [0,10], X0, korak, F, f1, C1, f2, C2, gradf1, gradf2);
    K=Y;
  % Izracun po RK4
  elseif(funkcija == 0)
    f = @(t,Y) [F(Y)];
    [t, Y] = rk4(f, [0,5], X0, korak, F, f1, C1, f2, C2, gradf1, gradf2);
    K=Y;
  endif
  % Narisemo obe funkciji ter presecisce
  [y,x,z] = ndgrid(linspace(-10,10,64));
  cla reset;
	f1 = @(x,y,z) e.^(-x.^2+1)+y.^2+z.^2;
        f2 = @(x,y,z) x.^2 + y.^2 + z.^2;

  f = f1(x,y,z);
  g = f2(x,y,z);
  isosurface(x,y,z,f,C1);
  hold on;
  isosurface(x,y,z,g,C2);
  axis equal;
  plot3(Y(1,:),Y(2,:),Y(3,:),"*r");
endfunction

function [t, Y] = euler(f, tspan, Y0, n, F,f1, C1, f2, C2, gradf1, gradf2)
  % [t, Y] = euler(f, tspan, Y0, n) poisce priblizno resitev DE 
  % dY/dt = f(t, Y) z zacetnim pogojem Y(t0) = Y0 
  % z Eulerjevo metodo z n koraki na intervalu tspan = [t0, tk].

  %h = (tspan(2) - tspan(1))/(n-1);
  h = 0.00001;
  t = linspace(tspan(1), tspan(2), n);
  Y = Y0;
  for k = 1:(n-1)
    
    y = Y(:,k) + feval(f, t(k), Y(:,k))*h;
    v = F(y); 
    JG = @(X)[gradf1(X)'; gradf2(X)'; v'];
    G = @(X) [f1(X)-C1; f2(X)-C2; v'*X-v'*y];
    [Y(:,k+1), stK ]= newton(G, JG, y, 10^-8, 10); % popravljanje priblizka
    
    %logika adaptivnega koraka
    if (stK <= 2)
      h = min(h * 2,100);
    else
      h = max(h / 2,0.0000001);
      if h != 0.0000001
      k = k - 1;
      endif
    endif  
    
  end
endfunction

function [t, Y] = rk4(f, tspan, Y0, n, F,f1, C1, f2, C2, gradf1, gradf2)
  % [t, Y] = rk4(f, tspan, Y0, n) poisce priblizno resitev DE 
  % dY/dt = f(t, Y) z zacetnim pogojem Y(t0) = Y0 
  % z Runge-Kutta metodo 4. reda z n koraki na intervalu 
  % tspan = [t0, tk].

  %h = (tspan(2) - tspan(1))/(n-1);
  h = 0.00001;
  t = linspace(tspan(1), tspan(2), n);
  Y = Y0;
  for k = 1:(n-1)
    k1 = h*feval(f, t(k), Y(:,k));
    k2 = h*feval(f, t(k) + h/2, Y(:,k) + k1/2);
    k3 = h*feval(f, t(k) + h/2, Y(:,k) + k2/2);
    k4 = h*feval(f, t(k) + h, Y(:,k) + k3);
    y = Y(:,k) + (k1 + 2*k2 + 2*k3 + k4)/6;
    v = F(y); 
    JG = @(X)[gradf1(X)'; gradf2(X)'; v'];
    G = @(X) [f1(X)-C1; f2(X)-C2; v'*X-v'*y];
    [Y(:,k+1),stK] = newton(G, JG, y, 10^-8, 10); % popravljanje priblizka
    
    %logika adaptivnega koraka
    if (stK <= 2)
      h = min(h * 2,100);
    else
      h = max(h / 2,0.0000001);
      if h != 0.0000001
      k = k - 1;
      endif
    endif  
    
  endfor
endfunction

function [x, k] = newton(G, JG, x0, tol, maxit)
  %[x, k] = newton(G, JG, x0, tol, maxit) poisce priblizek 
  %x za resitev enacbe F(x) = 0 z Newtonovo metodo.
  %(k je stevilo korakov, ki jih metoda porabi.)
  %G... funkcija
  %JG... Jacobijeva matrika G
  %x0... zacetni priblizek
  %tol... zahtevana natancnost
  %maxit... najvecje stevilo iteracij

  for k = 1:maxit
    %Izvedemo en korak Newtonove metode...
    %endif
    x = x0 - feval(JG, x0)\feval(G, x0);
    %... in testiramo, ce je metoda ze 'konvergirala'.
    if(norm(x - x0) < tol)
      break;
    end
    x0 = x;
  end
  %Izpisemo opozorilo v primeru, da zadnji priblizek ni znotraj tolerancnega obmocja.
  if(k == maxit)
    %disp("Warning: The method did not converge after maxit iterations.")
  end
endfunction

%!test
% [K]