close all
resynthesizeBool = 0;
if resynthesizeBool; clear all; resynthesizeBool = 1; end
saveFigBool = 0;

aRes = 13;
a = linspace(1,4,aRes)*5;
tStepHours = 12;
ctrlhorizonDays = 14;
T = ctrlhorizonDays*24/tStepHours;


J.norm = 1;
J.xtmp = 0;
J.utmp = 1;
J.uncertainEps = 0.0;
J.ifnCost = 0;
J.allDelay = 0/tStepHours;
p = 5;

halfLifeDays = 2;
k = 1-(.5)^(tStepHours/(24*halfLifeDays));

txf = @(x) sum(x);

allvcases = {'Ideal','Time-Varying','Time- and Host-Varying','Time- and Host-Varying with Population Stucture'};
vx = zeros(length(allvcases),aRes);
tx = zeros(length(allvcases),aRes);
timevx = zeros(aRes,10);
skiptx = zeros(aRes,10);

eclipseTime = 12;
eclipseTimeSteps = eclipseTime/tStepHours;

if resynthesizeBool; SynthesizeImmune; end

% Virulence
for vInd = 1:length(allvcases)
    SetPrefsV

    for ii = length(a):-1:1
        if ifnDelayedBool
            Phiu = allPhiui{ii,end};
            Phix = allPhixi{ii,end};
        else
            Phiu = allPhiui{ii,1};
            Phix = allPhixi{ii,1};
        end        
        xsy{ii,vInd} = Jo.umat(:,1:N:end).*Phiu(:,1:N:end);
        vx(vInd,ii) = norm(Jo.umat.*Phiu,J.norm);
%         vx(vInd,ii) = norm(xsy{ii,vInd},J.norm);

    end
end

%
vxe= vx(2:end,:);
vxrange = range(vxe(:));
vxmin = min(vxe(:));
gammafx = @(a) max(min((exp(-p*(a-vxmin)/vxrange)),1),0);
% Transmission

for ii = 1:length(a)
    A = allA{ii};
    for vInd = 1:length(allvcases)
        SetPrefsV
        if ifnDelayedBool
            Phiu = allPhiui{ii,end};
            Phix = allPhixi{ii,end};
        else
            Phiu = allPhiui{ii,1};
            Phix = allPhixi{ii,1};
        end
        if popStructBool
            rcva = a(ii);
        else
            rcva = median(a);
        end
        w = (qNorm([(1-skipRate)*vx(vInd,ii)*[1; zeros(N-1,1)] xsy{ii,vInd}],J.norm));
        shed = [(sum(Phix(:,1:N:end,:),1)) 0];
        fullts = shed.*gammafx(w);
        tx(vInd,ii) = txf(rcva*fullts);
        allShed{ii} = shed;
        allSym{ii} = w;
        allTx{ii} = fullts;

        if vInd == (length(allvcases))
            for ss = 1:10
                ws = qNorm(([(1-ss/10)*vx(vInd,ii)*[1; zeros(N-1,1)]  xsy{ii,vInd}]),J.norm);
                fullts = shed.*gammafx(ws);
                skiptx(ii,ss) = txf(rcva*fullts);
                
                if ss <= 5*24/tStepHours
                    timevx(ii,ss) = norm(Jo.umat.*allPhiui{ii,1+ss},J.norm);
                else
                    timevx(ii,ss) = nan;
                end
                
            end
        end
        enmOnly(ii) = norm(Jo.umat(:,1:N:end).*allPhiue{ii,1}(:,1:N:end),J.norm);
        enmTsImm(ii,:) = allPhiue{ii,2}(1,1:N:end);
        enmTsVir(ii,:) = allPhixe{ii,2}(1,1:N:end);
        ifnSupp(ii,:) = norm(Jo.umat(:,1:N:end).*allPhiui{ii,1+4*24/tStepHours}(:,1:N:end),J.norm);
        bestOfBoth(ii,:) = norm(Jo.umat(:,1:N:end).*allPhiui{ii,1}(:,1:N:end),J.norm);

        for tt = 1:(5*24/tStepHours)
            Ae = zeros(N);
            Ae(1,1) = (max(real(eig(A))));
            olOnly(ii,tt) = olfx(Ae,tt-1);
        end
    end
end
%
MakeVttPictures
