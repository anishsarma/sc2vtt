function outCurve =  GetOscEnvelope(inCurve)

outCurve=inCurve;

dips = find((diff(inCurve))<0)+1;
for dInd = 1:length(dips)
    tInd = dips(dInd);
    if tInd > 1 && tInd < length(outCurve)
        outCurve(tInd) = mean([outCurve(tInd-1) outCurve(tInd+1)]);
    end
end
% 
dips = find(diff(diff(outCurve))<0)+2;
for dInd = 1:length(dips)
    tInd = dips(dInd);
    if tInd > 1 && tInd < length(outCurve)
        outCurve(tInd) = mean([outCurve(tInd-1) outCurve(tInd+1)]);
    end
end

end