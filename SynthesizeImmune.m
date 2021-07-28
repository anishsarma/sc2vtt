for ii = 1:length(a)
    A = [1-k zeros(1,eclipseTimeSteps-1) k*a(ii);
        eye(eclipseTimeSteps) zeros(eclipseTimeSteps,1)];

    N = size(A,1);
    B = zeros(N);
    B(1,1) = 1;
    B(2,2) = 1;
    allA{ii} = A;
    for tt = 1:T
        J.x{tt} = J.xtmp * eye(N);
        J.u{tt} = zeros(N);
        J.u{tt}(1,1) = J.utmp;
        J.u{tt}(2,1) = J.utmp * J.ifnCost;
    end
    for tt = 1:11
    J.ifnDelay = tt-1; 
    J.allDelay = tt-1;

    [Phixi,Phiui,Jo] = SynthSLS(T,A,B,J);
    allPhixi{ii,tt} = Phixi;
    allPhiui{ii,tt} = Phiui;
    if any(any(isnan([Phixi,Phiui])))
        error('Nan')
    end
    end
    J.ifnDelay = T + 10;
    J.allDelay = 0;
    [Phixe,Phiue] = SynthSLS(T,A,B,J);
    allPhixe{ii,1} = Phixe;
    allPhiue{ii,1} = Phiue;
    if any(any(isnan([Phixe,Phiue])))
        error('Nan')
    end
    
    Ae = zeros(N);
    Ae(1,1) = (max(real(eig(A))));
    Be = eye(N);
    J.ifnDelay = T + 10;
    [Phixe,Phiue] = SynthSLS(T,Ae,Be,J);
    allPhixe{ii,2} = Phixe;
    allPhiue{ii,2} = Phiue;
    if any(any(isnan([Phixe,Phiue])))
        error('Nan')
    end


end