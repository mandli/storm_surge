      dlam = 1.e9;
      bb = pi*2*1.e8;
      dd = 2.e4;
      capf = 1.e0;
      rr = 0.02e0;
      ff = 1.e-13;
      gamma = capf*pi / (rr*bb);
      alfs = (dd/rr) * ff;
      capa = -alfs/2 + sqrt(alfs^2/4 + (pi/bb)^2);
      capb = -alfs/2 - sqrt(alfs^2/4 + (pi/bb)^2);
      pp = (1 - exp(capb*dlam)) / (exp(capa*dlam) - exp(capb*dlam));
      qq = 1 - pp;
      dgbp2 = gamma * (bb/pi)^2;

      dx = 7e6;
      [x,y] = meshgrid(0:dx:7e8, 0:dx:7e8);
      x0 = 1.5e8;
      y0 = 0;
      theta = .15;
      alf = cos(theta);
      beta = sin(theta);
      xbar = alf*(x-x0) + beta*(y-y0);
      ybar = -beta*(x-x0) + alf*(y-y0);
      stream = -dgbp2 * sin(pi*ybar/bb) .* ...
             (pp*exp(capa*xbar) + qq*exp(capb*xbar) - 1);
      stream(find(xbar<0)) = 0;	
      stream(find(ybar>pi*2e8)) = 0;	
      cline = -4e8:1e8:2e9;
      contour(x,y,stream,cline)
      axis square
