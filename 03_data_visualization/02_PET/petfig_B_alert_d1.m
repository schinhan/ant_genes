%% makes the pet-plot part A
% requirements:
% ENIGMA toolbox for plotting:
 system(['git clone https://github.com/MICA-MNI/ENIGMA.git'])
% altmany's export figure tool
system(['git clone https://github.com/altmany/export_fig.git'])

%% prepare
% load required info on surface 
load visualization/pet/plotstuff.mat

% load maps to plot (all in correct cato-order)
load visualization/pet/stuff4plot.mat

% prerpare parcellation for surface
deg_fsa5(:,1) = parcel_to_surface([0;zscore(ant_recep(:,targets(2)))], 'laus_250_fsa5');
deg_fsa5(:,2) = parcel_to_surface([0;zscore(ant_recep(:,targets(2)))], 'laus_250_fsa5');

for i=1:3
    antpl(i,:)=parcel_to_surface([0;ant_cons(:,i)], 'laus_250_fsa5');
end

%% then plot
figure;
%% plot ANT maps
subplot('position',[.28 .81 .16 .20])
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
    double(antpl(1,vl)),'EdgeColor','none');
    view(-90,0)
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;
    
subplot('position',[.45 .81 .16 .20])
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
    double(antpl(1,vl)),'EdgeColor','none');
    view(90,0)
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;
   
subplot('position',[.62 .81 .16 .20])
    trisurf(surf.tri(tr,:)-length(vr),surf.coord(1,vr),surf.coord(2,vr),surf.coord(3,vr),...
    double(antpl(1,vr)),'EdgeColor','none');
    view(-90,0); 
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;
subplot('position',[.79 .81 .16 .20])
    trisurf(surf.tri(tr,:)-length(vr),surf.coord(1,vr),surf.coord(2,vr),surf.coord(3,vr),...
    double(antpl(1,vr)),'EdgeColor','none');
    view(90,0); 
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;

%% plot receptor maps

subplot('position',[.02 .62 .16 .20])
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
    double(deg_fsa5(vl,1)),'EdgeColor','none');
    view(-90,0)
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;
    
subplot('position',[.02 .42 .16 .20])
    trisurf(surf.tri(tl,:),surf.coord(1,vl),surf.coord(2,vl),surf.coord(3,vl),...
    double(deg_fsa5(vl,1)),'EdgeColor','none');
    view(90,0)
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;
    
subplot('position',[.02 .22 .16 .20])
    trisurf(surf.tri(tr,:)-length(vr),surf.coord(1,vr),surf.coord(2,vr),surf.coord(3,vr),...
    double(deg_fsa5(vr,1)),'EdgeColor','none');
    view(-90,0); 
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;
subplot('position',[.02 .02 .16 .20])
    trisurf(surf.tri(tr,:)-length(vr),surf.coord(1,vr),surf.coord(2,vr),surf.coord(3,vr),...
    double(deg_fsa5(vr,1)),'EdgeColor','none');
    view(90,0); 
    daspect([1 1 1]); axis tight; camlight; axis vis3d off;
    lighting phong; material dull; shading flat;

   
%% plot scatter 
subplot('position',[.31 .13 .63 .65])
    scatter(ant_cons(:,1),ant_recep(:,targets(2)),40,'k','filled');l=lsline; l.Color='r'; l.LineWidth=2;
    xlabel('ANT Alerting'); ylabel('D1 availability')
%% adjust colormaps
kids=get(gcf,'Children');
colormap('inferno');
for i=1:length(kids)
    klim(i,:)=kids(i).CLim;kpos(i,:)=kids(i).Position;
end

for i=2:5
    kids(i).CLim=[min(klim(kpos(:,1)==.02,1)) max(klim(kpos(:,1)==.02,2))];
   
end
for i=6:9
    kids(i).CLim=[min(klim(kpos(:,2)==.81,1)) max(klim(kpos(:,2)==.81,2))];
end

%% export            
set(gcf,'color','w');
export_fig('visualization/pet/petfig_B_alert_d1','-pdf') 

