%%

newColors = [31 135 198;
   239 139 38;
     231 32 76;]/255;
 
pltThese = [1 ceil(aRes/2) aRes];
txScale = 2/min(max(tx([2  3 4],:),[],2));
lw = 4;
figure;

jj = 1;
for ii = 1:size(olOnly,1)
    if any(ii==pltThese)
        pcol = newColors(jj,:);
        jj= jj + 1;
        lw2 = lw;
    elseif mod(ii,2) == 1
        pcol = [.8 .8 .8];
        lw2 = lw / 2;
    else
        continue
    end
    tax=((1:length(olOnly(ii,:))))*tStepHours/24;
    semilogy(tax,olOnly(ii,:),'Color',pcol,'LineWidth',lw2);hold on;
    
end
xlabel('Time (days)')
ylabel('Viral Load without Control')
PictureFormat
axis([-Inf Inf -Inf Inf])
xlabel('Time (days)')
ylabel('Open-Loop Shedding')
SaveDownFigs(gcf,'OpenLoop',saveFigBool)


figure;
set(gcf,'defaultAxesColorOrder',newColors)
for jj = 1:length(pltThese)
     ii = pltThese(jj);
    shed = allShed{ii};
    w = allSym{ii};
    fullts = allTx{ii};
    tax= (1:length(shed))/24*tStepHours;
%     plot(tax,GetOscEnvelope(shed),'LineWidth',lw)
    plot(tax,(shed),'LineWidth',lw)
    hold on
    PictureFormat
    axis([-Inf Inf -Inf Inf])
end
xlabel('Time (days)')
ylabel('Shedding')

SaveDownFigs(gcf,'ShedSymTx1',saveFigBool)
figure;set(gcf,'defaultAxesColorOrder',newColors)

for jj = 1:length(pltThese)
     ii = pltThese(jj);
    shed = allShed{ii};
    w = allSym{ii};
    fullts = allTx{ii};
    tax= (1:length(w))/24*tStepHours;
    plot(tax,1-gammafx(w),'LineWidth',lw)
    hold on
    PictureFormat
    axis([-Inf Inf 0 1])

end
xlabel('Time (days)')
ylabel('Avoidance (1-\gamma)')

SaveDownFigs(gcf,'ShedSymTx2',saveFigBool)
figure;set(gcf,'defaultAxesColorOrder',newColors)

for jj = 1:length(pltThese)
     ii = pltThese(jj);
    shed = allShed{ii};
    w = allSym{ii};
    fullts = allTx{ii};
    tax= (1:length(fullts))/24*tStepHours;
    plot(tax,cumsum(a(ii)*fullts)*txScale,'LineWidth',lw)
    hold on
    PictureFormat
    axis([-Inf Inf -Inf Inf])

end
xlabel('Time (days)')
ylabel('Cumulative Transmission')

SaveDownFigs(gcf,'ShedSymTx3',saveFigBool)


figure;set(gcf,'defaultAxesColorOrder',newColors)

pctCell = 1:1:100;
maxFc = @(x) 1+(100-x)./x;
onepFc = @(x) (min(x+1,100))./x;
plot(pctCell,maxFc(pctCell),'k','LineWidth',lw)
hold on;
plot(pctCell,onepFc(pctCell),'k--','LineWidth',1)
axis([-Inf Inf 0 12])
hold on
xlabel('Median SCP')
ylabel('Fold Change')
PictureFormat
SaveDownFigs(gcf,'Sparsity',saveFigBool)

figure;set(gcf,'defaultAxesColorOrder',newColors)

for ii = 1:3
    tax= ((1:size(enmTsImm,2)))/24*tStepHours;
    plot(tax,(enmTsImm(pltThese(ii),:)),'LineWidth',lw); hold on;
end
xlabel('Time (days)')
ylabel('Immune Effort')
axis([-Inf Inf -Inf Inf])
PictureFormat
SaveDownFigs(gcf,'MoreEffortSameOutcome1',saveFigBool)

figure;set(gcf,'defaultAxesColorOrder',newColors)

for ii = 1:3
    tax= ((1:size(enmTsImm,2)))/24*tStepHours;
    plot(tax,enmTsVir(pltThese(ii),:),'LineWidth',lw); hold on;
end
xlabel('Time (days)')
ylabel('Controlled Viral Load')

axis([-Inf Inf -Inf Inf])
PictureFormat
SaveDownFigs(gcf,'MoreEffortSameOutcome2',saveFigBool)

figure; set(gcf,'defaultAxesColorOrder',newColors)

plot((1:10)/10,skiptx(pltThese,:)'*txScale,'LineWidth',lw);
xlabel('Escape Rate')
ylabel('Effective Transmission')

axis([-Inf Inf -Inf Inf])
PictureFormat
SaveDownFigs(gcf,'AsymEscape',saveFigBool)

figure; set(gcf,'defaultAxesColorOrder',newColors)

plot(linspace(1*tStepHours/24,10*tStepHours/24,10),timevx(pltThese,:)','LineWidth',lw);
axis([-Inf Inf -Inf Inf])

ylabel('Virulence (a.u.)')
xlabel('Time to Interferon (days)')
PictureFormat
SaveDownFigs(gcf,'InterferonTime',saveFigBool)

figure;
for vInd = 1:length(allvcases)
    
    plot(a,bestOfBoth,'-','Color','#22B573','LineWidth',lw); hold on;

    plot(a,enmOnly(1,:),'Color','#A21D31','LineWidth',lw); hold on;
    plot(a,ifnSupp,'-','Color','#FF99FF','LineWidth',lw); hold on;
    xlabel('Host-Viral \alpha')
    ylabel('Virulence (a.u.)')
    axis([-Inf Inf -Inf Inf])
end

PictureFormat
SaveDownFigs(gcf,'VirulenceAlpha',saveFigBool)

newcolors = {'#22B573','#EF8B26','#9856A1','#9856A1'};

figure;
for vInd = 1:length(allvcases)
    if vInd == length(allvcases)
        ls = '--';
    else
        ls = '-';
    end
    if vInd == 1 
        hp = plot(vx(vInd,1),tx(vInd,1).*txScale,'o','Color',newcolors{vInd}, 'LineWidth',lw,'MarkerFaceColor','b','MarkerSize',16); hold on
        set(hp,'MarkerFaceColor',get(hp,'Color'))
    else
        hp = plot(vx(vInd,:),tx(vInd,:).*txScale,ls,'Color',newcolors{vInd},'LineWidth',lw); hold on
    end
    ylabel('Transmission');
    xlabel('Virulence (a.u.)')
    axis([-Inf Inf -Inf Inf])
end
PictureFormat
SaveDownFigs(gcf,'VirulenceTransmission',saveFigBool)
