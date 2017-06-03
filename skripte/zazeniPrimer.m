function K = zazeniPrimer(primer,korak, funkcija,adaptivno)
  %primer     ...oznacuje steviko primera
  % korak    ... Korak pri Eulerjevi oz RK4 metodi
  % funkcija ... Ce 1: Euler
  %              Ce 0: RK4
  % adaptivno ... Ce 1: adaptivno
  %              Ce 0: fiksno
  C1 = 4;
  C2 = 1;
  f1 = NaN;
  f2= NaN;
  gradf1= NaN;
  gradf1= NaN;
  pf1= NaN;
  pf2= NaN;
  
  
  X0=[1;1;sqrt(2)];
  switch (primer)
  case 1
% ************************************* %
% DELUJOC PRIMER 1 sfera in valj
% ************************************* %
  f1 = @(a) a(1).^2 + a(2).^2 + a(3).^2;
  gradf1 = @(a) [2*a(1);2*a(2);2*a(3)];

	f2 = @(a) a(1).^2 + a(2).^2;
  gradf2 = @(a) [2*a(1);2*a(2);0];

  % IZRIS
	pf1 = @(x,y,z) x.^2 + y.^2 + z.^2;
  pf2 = @(x,y,z) x.^2 +y.^2;

  case 2
    % ************************************* %
% DELUJOC PRIMER 2, sfera, ravnina
% ************************************* %

	f1 = @(a) a(1).^2 + a(2).^2 + a(3).^2;
  gradf1 =@(a) [2*a(1);2*a(2);2*a(3)];
	
	f2 = @(a) a(1)*3 + a(2)*2 + a(3);
  gradf2 = @(a) [3; 2; 1];

	% IZRIS
	pf1 = @(x,y,z) x.^2 + y.^2 + z.^2;
	pf2 = @(x,y,z) 3*x + 2*y + z;

  case 3
      % ************************************* %
  % DELUJOC PRIMER 3, sfera in 2: iz pdf
  % ************************************* %

  f1 = @(a) a(1).^2 + a(2).^2 + a(3).^2;
  gradf1 =@(a) [2*a(1);2*a(2);2*a(3)];

	f2 = @(a) a(2).^4 + log(a(1).^2 + 1).*a(3).^2 - 4;
  gradf2 = @(a) [(2*a(1).*a(3).^2)/(a(1).^2+1); 4*a(2).^3; 2*a(3).*log(a(1).^2 + 1)];

	% IZRIS
	pf1 = @(x,y,z) x.^2 + y.^2 + z.^2;
	pf2 = @(x,y,z) y.^4 + log(x.^2 + 1).*z.^2 - 4;

  case 4
      % ************************************* %
  % DELUJOC PRIMER 4, 1: in 2: iz pdf
  % ************************************* %
    st korakov = 300;
    tspan = [0, 100]

    f1 = @(a) a(1).^2 + (cos(a(2)).*(a(3).^2)) - 1;
    gradf1 = @(a) [2*a(1); -a(3).^2.*sin(a(2)); 2*a(3).*cos(a(2))];

    f2 = @(a) a(2).^4 + log(a(1).^2 + 1).*a(3).^2 - 4;
    gradf2 = @(a) [(2*a(1).*a(3).^2)/(a(1).^2+1); 4*a(2).^3; 2*a(3).*log(a(1).^2 + 1)];

    % IZRIS
    pf1 = @(x,y,z) x.^2 + cos(y).*z.^2 - 1;
    pf2 = @(x,y,z) y.^4 + log(x.^2 + 1).*z.^2 - 4;
    C2 = 4*pi*cos(45);

  case 5
   
    %Tukaj sem vzela funkcije 3 in 5 v testi (C1=3,C2=10)
    C1 = 3;
    C2 = 10;

   f1=@(a) e.^(-a(1).^2+1) + a(2).^2+a(3).^ 2;
   f2=@(a) e.^(a(1).*a(2).*a(3))+a(2).^2+a(3).^2;
   
   gradf1 =@(a) [-2*a(1)*e.^(-a(1).^2+1);2*a(2);2*a(3)];
   gradf2 =@(a) [a(2)*a(3)*e.^(a(1).*a(2).*a(3));a(1)*a(3)*e.^(a(1).*a(2).*a(3))+2*a(2);a(1)*a(2)*e.^(a(1).* a(2).*a(3))+2*a(3)];

  % Ko izrisem:
	pf1=@(x,y,z) e.^(-x.^2+1)+y.^2+z.^2;
   	pf2=@(x,y,z) e.^(x.*y.*z)+y.^2+z.^2;

    
  case 6
    %Tukaj sem vzela funkcije 3 in sfera (C1=3,C2=4)
    C1 = 3;
    C2 = 4;
   f1 = @(a) e.^(-a(1).^2+1) + a(2).^2+a(3).^ 2;
   f2 = @(a) a(1).^2 + a(2).^2 + a(3).^2;

   gradf1 =@(a) [-2*a(1)*e.^(-a(1).^2+1);2*a(2);2*a(3)];
   gradf2 =@(a) [2*a(1);2*a(2);2*a(3)];

   %Ko izrisem:
	pf1 = @(x,y,z) e.^(-x.^2+1)+y.^2+z.^2;
   pf2 = @(x,y,z) x.^2 + y.^2 + z.^2;
  case 7
    %Tukaj sem vzela funkcije 5 in sfera (C1=10,C2=4)
    C1 = 10;
    C2 = 4;
   f1=@(a) e.^(a(1).*a(2).*a(3))+a(2).^2+a(3).^2;
   f2 = @(a) a(1).^2 + a(2).^2 + a(3).^2;

   gradf1 =@(a) [a(2)*a(3)*e.^(a(1).*a(2).*a(3));a(1)*a(3)*e.^(a(1).*a(2).*a(3))+2*a(2);a(1)*a(2)*e.^(a(1).* a(2).*a(3))+2*a(3)];
   gradf2 =@(a) [2*a(1);2*a(2);2*a(3)];

  % Ko izrisem:
	  pf1=@(x,y,z) e.^(x.*y.*z)+y.^2+z.^2;
   pf2 = @(x,y,z) x.^2 + y.^2 + z.^2;
  case 8
%*Tukaj sem vzela funkcije 3 in valj(C1=3,C2=1)
    C1 = 3;
    C2 = 1;
   f1 = @(a) e.^(-a(1).^2+1) + a(2).^2+a(3).^ 2;
   f2 = @(a) a(1).^2 + a(2).^2;

   gradf1 =@(a) [-2*a(1)*e.^(-a(1).^2+1);2*a(2);2*a(3)];
   gradf2 =@(a) [2*a(1);2*a(2);0];

   %Ko izrisem:
	pf1 = @(x,y,z) e.^(-x.^2+1)+y.^2+z.^2;
	pf2 = @(x,y,z) x.^2 +y.^2;   
  case 9
    %*Tukaj sem vzela funkcije 5 in valj(C1=10,C2=1)
    C1 = 10;
    C2 = 1;
   f1=@(a) e.^(a(1).*a(2).*a(3))+a(2).^2+a(3).^2;
   f2 = @(a) a(1).^2 + a(2).^2;

   gradf1 =@(a) [a(2)*a(3)*e.^(a(1).*a(2).*a(3));a(1)*a(3)*e.^(a(1).*a(2).*a(3))+2*a(2);a(1)*a(2)*e.^(a(1).* a(2).*a(3))+2*a(3)];
   gradf2 =@(a) [2*a(1);2*a(2);0];

  % Ko izrisem:
	pf1=@(x,y,z) e.^(x.*y.*z)+y.^2+z.^2;
	pf2 = @(x,y,z) x.^2 +y.^2;
endswitch
  
  
  if adaptivno == 1
  K = presekDvehPloskevAdaptTest(f1, C1, f2, C2, gradf1, gradf2, X0, korak, funkcija,pf1,pf2);
  else
  K = presekDvehPloskevTest(f1, C1, f2, C2, gradf1, gradf2, X0, korak, funkcija,pf1,pf2);
  endif
  
  
endfunction