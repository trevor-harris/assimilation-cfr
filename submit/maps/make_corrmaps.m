
%===========================
% CORRELATION MAPS OF PHYDA
%===========================

addpath('../')

% ALL PROXY TYPES, LMR DATABASE
%efl='cesm_lme010_r12000_p998_state_t2qjozn_avg_JunAug_prxydtst_5_prxtp_tca_2896_swtchbld100_05-Jan-2018_14:49:45.mat';
efl='cesm_lme010_r12000_p998_state_t2qjozn_avg_AprMar_prxydtst_5_prxtp_tca_2896_swtchbld100_05-Jan-2018_14:33:07.mat';
% INDEX ONLY RECONSTRUCTION WITH FULL ENSEMBLES
%efl='cesm_lme010_r12000_p998_state_2ozn_avg_AprMar_prxydtst_5_prxtp_tca_2605_swtchbld100_16-Nov-2017_15:57:59.mat';


%pth='/d2/nsteiger/output-da/hydroclimate/'; load([pth,efl])
pth='/d2/nsteiger/output-da/hydroclimate/subrecons_5Jan18/';load([pth,efl])


%varpth=['/d6/haibo/lme/cam5/BLMTRC5CN/ALLFORCEINGS/010'];
%sftlf=permute(ncread([varpth '/QSOIL/b.e11.BLMTRC5CN.f19_g16.010.clm2.h0.QSOIL.085001-184912.nc'],'landmask'),[2 1]);
%lat=ncread([varpth '/QSOIL/b.e11.BLMTRC5CN.f19_g16.010.clm2.h0.QSOIL.085001-184912.nc'],'lat');
%lon=ncread([varpth '/QSOIL/b.e11.BLMTRC5CN.f19_g16.010.clm2.h0.QSOIL.085001-184912.nc'],'lon');

ens='010';
varpth=['/d6/haibo/lme/cam5/BLMTRC5CN/ALLFORCEINGS/' ens];
lat=ncread([varpth '/TREFHT/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TREFHT.085001-184912.nc'],'lat');
lon=ncread([varpth '/TREFHT/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TREFHT.085001-184912.nc'],'lon');

% Surface temperature
Xv_mon=ncread([varpth '/TREFHT/b.e11.BLMTRC5CN.f19_g16.' ens '.cam.h0.TREFHT.085001-184912.nc'],'TREFHT');
Xv_mon=permute(Xv_mon,[2 1 3]);
cesm_t2m_a2m=mon2ann(Xv_mon,mon_avg_o,mon_avg_f);
   		
[cesm_nino34] = nino(cesm_t2m_a2m,lat,lon,3,'n');

% PDSI
load(['/d1/nsteiger/climate-data/cesm-lme/cesm_lme_' ens '_pdsi_output_AWC_c_7_u2_gcm.mat'],'pdsi_f')
cesm_pdsi_jja=mon2ann(pdsi_f,6,8);
cesm_pdsi_djf=mon2ann(pdsi_f,12,2);


%-------------------------
% LOAD SPATIAL FIELDS
%-------------------------

% % Temperature
% k=find(state_tp=='t');
% Xt2m=Xa_m(xmeta{k}.id_X(1):xmeta{k}.id_X(2),:);
% %Xt2m=reshape(Xvar,length(xmeta{k}.lat),length(xmeta{k}.lon),reconYrs);
% 
% 
% clearvars -except xa_amo xb_itcz_sa xb_itcz_atl xb_itcz_pac xa_itcz_sa xa_itcz_atl xa_itcz_pac Xt2m xa_nino xb_nino
% 
% 
% % JJA PDSI
% efl='cesm_lme010_r12000_p998_state_t2qjozn_avg_JunAug_prxydtst_5_prxtp_tca_2896_swtchbld100_05-Jan-2018_14:49:45.mat';
% pth='/d2/nsteiger/output-da/hydroclimate/subrecons_5Jan18/';load([pth,efl])
% %pth='/d2/nsteiger/output-da/hydroclimate/'; load([pth,efl])


% What years to look at?
analyrs1=855:1848;
%analyrs2=1601:1925;
%analyrs3=800:2000;
%analyrs4=800:1925;
[~,~,ir1] = intersect(analyrs1,850:1848); % indices of recon segment
%[~,~,ir2] = intersect(analyrs2,r_o:r_f); % indices of recon segment
%[~,~,ir3] = intersect(analyrs3,r_o:r_f); % indices of recon segment
%[~,~,ir4] = intersect(analyrs4,r_o:r_f); % indices of recon segment
% monthly indices
ovrlprng=rng_intr([r_o,r_f],[analyrs1(1),analyrs1(end)]);
mon_r=r_o:(1/12):(r_f+11/12);
irm1_a=find(abs(ovrlprng(1)-mon_r)<1e-6);
irm1_b=find(abs(ovrlprng(2)-1/12-mon_r)<1e-6)+12;

% ovrlprng=rng_intr([r_o,r_f],[analyrs2(1),analyrs2(end)]);
% mon_r=r_o:(1/12):(r_f+11/12);
% irm2_a=find(abs(ovrlprng(1)-mon_r)<1e-6);
% irm2_b=find(abs(ovrlprng(2)-1/12-mon_r)<1e-6)+12;
% 
% ovrlprng=rng_intr([r_o,r_f],[analyrs3(1),analyrs3(end)]);
% mon_r=r_o:(1/12):(r_f+11/12);
% irm3_a=find(abs(ovrlprng(1)-mon_r)<1e-6);
% irm3_b=find(abs(ovrlprng(2)-1/12-mon_r)<1e-6)+12;


% % PDSI
% k=find(state_tp=='q');
% Xvar=Xa_m(xmeta{k}.id_X(1):xmeta{k}.id_X(2),:);
% Xpdsi=NaN(length(xmeta{k}.lat)*length(xmeta{k}.lon),reconYrs);
% Xpdsi(xmeta{k}.lndidx,:)=Xvar;
% Xpdsi_jja=Xpdsi;
% Xpdsi=reshape(Xpdsi,length(xmeta{k}.lat),length(xmeta{k}.lon),reconYrs);

% lat=xmeta{k}.lat;
% lon=xmeta{k}.lon;



% % PULL OUT SPATIAL MEAN OVER AMERICAN SOUTHWEST
% %load('/d1/nsteiger/climate-data/cesm-lme/Nrth_A_Swthwst_mask.mat')
% varpth=['/d6/haibo/lme/cam5/BLMTRC5CN/ALLFORCEINGS/010'];
% sftlf=permute(ncread([varpth '/QSOIL/b.e11.BLMTRC5CN.f19_g16.010.clm2.h0.QSOIL.085001-184912.nc'],'landmask'),[2 1]);
% %lat=ncread([varpth '/QSOIL/b.e11.BLMTRC5CN.f19_g16.010.clm2.h0.QSOIL.085001-184912.nc'],'lat');
% %lon=ncread([varpth '/QSOIL/b.e11.BLMTRC5CN.f19_g16.010.clm2.h0.QSOIL.085001-184912.nc'],'lon');
% % NASW = box of 235:255 E (125:105 W); 31:42 N (just the US, not Mexico)
% % Find nearest lat/lon equal to or outside bounds
% lonext=[235 255];latext=[31 42];
% ln_a=find((abs(lon-lonext(1))-min(abs(lon-lonext(1))))<1e-5);
% ln_b=find((abs(lon-lonext(2))-min(abs(lon-lonext(2))))<1e-5);
% lt_a=find((abs(lat-latext(1))-min(abs(lat-latext(1))))<1e-5);
% lt_b=find((abs(lat-latext(2))-min(abs(lat-latext(2))))<1e-5);
% nasw=zeros(size(sftlf));
% nasw(lt_a:lt_b,ln_a:ln_b)=sftlf(lt_a:lt_b,ln_a:ln_b);
% 
% nasw(nasw==0)=NaN;
% Xv_msk=bsxfun(@times,Xpdsi,nasw);
% A=cosd(repmat(xmeta{k}.lat,[1 length(xmeta{k}.lon)]));
% xa_di_jja=wmean_a(Xv_msk,A); % RECONSTRUCTION DROUGHT INDEX



% % DJF PDSI
% efl='cesm_lme010_r12000_p998_state_t2qjozn_avg_DecFeb_prxydtst_5_prxtp_tca_2896_swtchbld100_05-Jan-2018_14:40:23.mat';
% pth='/d2/nsteiger/output-da/hydroclimate/subrecons_5Jan18/';load([pth,efl])
% 
% % PDSI
% k=find(state_tp=='q');
% Xvar=Xa_m(xmeta{k}.id_X(1):xmeta{k}.id_X(2),:);
% Xpdsi=NaN(length(xmeta{k}.lat)*length(xmeta{k}.lon),reconYrs);
% Xpdsi(xmeta{k}.lndidx,:)=Xvar;
% Xpdsi_djf=Xpdsi;
% Xpdsi=reshape(Xpdsi,length(xmeta{k}.lat),length(xmeta{k}.lon),reconYrs);
% 
% Xv_msk=bsxfun(@times,Xpdsi,nasw);
% xa_di_djf=wmean_a(Xv_msk,A); % RECONSTRUCTION DROUGHT INDEX
% 
% 
% % Compute nino index from temperature field (so ages are consistent with other fields)
% [xa_nino34] = nino(reshape(Xt2m,length(lat),length(lon),reconYrs),lat,lon,3,'n');


% Compute the correlations
tcr=zeros(length(lat)*length(lon),1);
%ptcr=zeros(length(lat)*length(lon),1);
%ptcr_d=zeros(length(lat)*length(lon),1);
pcr_jja=zeros(length(lat)*length(lon),1);
pcr_djf=zeros(length(lat)*length(lon),1);
cesm_t2m_a2m=reshape(cesm_t2m_a2m,length(lat)*length(lon),size(cesm_t2m_a2m,3));
cesm_pdsi_jja=reshape(cesm_pdsi_jja,length(lat)*length(lon),size(cesm_pdsi_jja,3));
cesm_pdsi_djf=reshape(cesm_pdsi_djf,length(lat)*length(lon),size(cesm_pdsi_djf,3));

for i=1:length(lat)*length(lon)
   %tcr(i)=corr(xa_nino34(ir4),Xt2m(i,ir4)');
   %ptcr(i)=corr(xa_di_jja(ir4),Xt2m(i,ir4)');
   % %ptcr_d(i)=corr(ann2xt(xa_di_jja(ir4),10),ann2xt(Xt2m(i,ir4)',10));
   %pcr_jja(i)=corr(xa_nino34(ir4),Xpdsi_jja(i,ir4)');
   %pcr_djf(i)=corr(xa_nino34(ir4),Xpdsi_djf(i,ir4)');
   
   tcr(i)=corr(cesm_nino34(ir1),cesm_t2m_a2m(i,ir1)');
   pcr_jja(i)=corr(cesm_nino34(ir1),cesm_pdsi_jja(i,ir1)');
   pcr_djf(i)=corr(cesm_nino34(ir1),cesm_pdsi_djf(i,ir1)');
end


%---------------------------
%  SAVE PATTERNS TO NETCDF
%---------------------------

savenc='y';
if savenc=='y';

   inputf='cesm_latlon.nc';

   % define the output variable
   %outputf=['corrMaps_PHYDA_nino_t2m_pdsi.nc'];
   outputf=['corrMaps_CESM10_nino_t2m_pdsi.nc'];
   %outputf='test_output.nc';inputf='ccsm4_latlon.nc';
   % grab lat/lon from model, put it into the output file
   eval(['!ncks -A -v lat,lon ',inputf,' ',outputf]);

   cp_tcr=permute(reshape(tcr,length(lat),length(lon)),[2 1]);
   cp_tcr(cp_tcr==NaN)=1e20;
   nm1=['corr_nino_t2m'];

   cp_pjja=permute(reshape(pcr_jja,length(lat),length(lon)),[2 1]);
   cp_pjja(cp_pjja==NaN)=1e20;
   nm2=['corr_nino_pdsiJJA'];

   cp_pdjf=permute(reshape(pcr_djf,length(lat),length(lon)),[2 1]);
   cp_pdjf(cp_pdjf==NaN)=1e20;
   nm3=['corr_nino_pdsiDJF'];

   nccreate(outputf,nm1,'Dimensions',{'lon',length(lon),'lat',length(lat)},'FillValue',1e20);
   ncwrite(outputf,nm1,cp_tcr);
   ncwriteatt(outputf,nm1,'missing_value',1e20);

   nccreate(outputf,nm2,'Dimensions',{'lon',length(lon),'lat',length(lat)},'FillValue',1e20);
   ncwrite(outputf,nm2,cp_pjja);
   ncwriteatt(outputf,nm2,'missing_value',1e20);

   nccreate(outputf,nm3,'Dimensions',{'lon',length(lon),'lat',length(lat)},'FillValue',1e20);
   ncwrite(outputf,nm3,cp_pdjf);
   ncwriteatt(outputf,nm3,'missing_value',1e20);


end





plots='n';
if plots=='y'

   s=load('coast');mlat=s.lat;mlon=s.long;
   load('./colormaps/NCV_blue_red.mat')

   varZ=reshape(tcr,length(lat),length(lon));

   figure
   hold on
   %h=axesm('MapProjection','wagner4','MapLatLimit',[(min(xlat)) (max(xlat))],'MapLonLimit',[min(xlon) max(xlon)]);
   h=axesm('MapProjection','eqdcylin','MapLatLimit',[(min(lat)) (max(lat))],'MapLonLimit',[min(lon) max(lon)]);
   %h=axesm('MapProjection','eqdcylin','MapLatLimit',[-60 max(lat)]);
   pcolorm(lat,lon,varZ); % add additional value to avoid line in pcolor plots
   plotm(mlat,mlon,'color',[0.5 0.5 0.5])
   colormap(cmap) % colormap from 'WhiteYellowOrangeRed'
   caxis([-1 1])
   colorbar
   %set(gca,'fontsize',18)
   title('corr(nino3.4, T2m)')
   %set(get(colorbar('location','southoutside'),'xlabel'),'string',sn{j},'fontsize',14)
   %set(get(colorbar('location','southoutside'),'xlabel'),'string','Correlation','fontsize',18)
   %set(cbh,'YTick',[-1:0.2:1])
   %delete( cbh )
   hold off

   %print(['figs/corrmap_nino_t2m_phyda.png'],'-dpng','-r300');

   varZ=reshape(pcr_jja,length(lat),length(lon));

   figure
   hold on
   %h=axesm('MapProjection','wagner4','MapLatLimit',[(min(xlat)) (max(xlat))],'MapLonLimit',[min(xlon) max(xlon)]);
   h=axesm('MapProjection','eqdcylin','MapLatLimit',[-60 max(lat)],'MapLonLimit',[min(lon) max(lon)]);
   pcolorm(lat,lon,varZ); % add additional value to avoid line in pcolor plots
   plotm(mlat,mlon,'color',[0.5 0.5 0.5])
   colormap(cmap) % colormap from 'WhiteYellowOrangeRed'
   caxis([-1 1])
   colorbar
   %set(gca,'fontsize',18)
   title('corr(nino3.4, JJA PDSI)')
   %set(get(colorbar('location','southoutside'),'xlabel'),'string',sn{j},'fontsize',14)
   %set(get(colorbar('location','southoutside'),'xlabel'),'string','Correlation','fontsize',18)
   %set(cbh,'YTick',[-1:0.2:1])
   %delete( cbh )
   hold off

   %print(['figs/corrmap_nino_jjapdsi_phyda.png'],'-dpng','-r300');


   varZ=reshape(pcr_djf,length(lat),length(lon));

   figure
   hold on
   %h=axesm('MapProjection','wagner4','MapLatLimit',[(min(xlat)) (max(xlat))],'MapLonLimit',[min(xlon) max(xlon)]);
   h=axesm('MapProjection','eqdcylin','MapLatLimit',[-60 max(lat)],'MapLonLimit',[min(lon) max(lon)]);
   pcolorm(lat,lon,varZ); % add additional value to avoid line in pcolor plots
   plotm(mlat,mlon,'color',[0.5 0.5 0.5])
   colormap(cmap) % colormap from 'WhiteYellowOrangeRed'
   caxis([-1 1])
   colorbar
   %set(gca,'fontsize',18)
   title('corr(nino3.4, DJF PDSI)')
   %set(get(colorbar('location','southoutside'),'xlabel'),'string',sn{j},'fontsize',14)
   %set(get(colorbar('location','southoutside'),'xlabel'),'string','Correlation','fontsize',18)
   %set(cbh,'YTick',[-1:0.2:1])
   %delete( cbh )
   hold off

   %print(['figs/corrmap_nino_djfpdsi_phyda.png'],'-dpng','-r300');


   varZ=reshape(ptcr,length(lat),length(lon));

   figure
   hold on
   %h=axesm('MapProjection','wagner4','MapLatLimit',[(min(xlat)) (max(xlat))],'MapLonLimit',[min(xlon) max(xlon)]);
   h=axesm('MapProjection','eqdcylin','MapLatLimit',[(min(lat)) (max(lat))],'MapLonLimit',[min(lon) max(lon)]);
   pcolorm(lat,lon,varZ); % add additional value to avoid line in pcolor plots
   plotm(mlat,mlon,'color',[0.5 0.5 0.5])
   colormap(cmap) % colormap from 'WhiteYellowOrangeRed'
   caxis([-1 1])
   colorbar
   %set(gca,'fontsize',18)
   title('corr(NASW PDSI, T2m)')
   %set(get(colorbar('location','southoutside'),'xlabel'),'string',sn{j},'fontsize',14)
   %set(get(colorbar('location','southoutside'),'xlabel'),'string','Correlation','fontsize',18)
   %set(cbh,'YTick',[-1:0.2:1])
   %delete( cbh )
   hold off

   %print(['figs/corrmap_naswpdsi_t2m_phyda.png'],'-dpng','-r300');



   varZ=reshape(ptcr_d,length(lat),length(lon));

   figure
   hold on
   %h=axesm('MapProjection','wagner4','MapLatLimit',[(min(xlat)) (max(xlat))],'MapLonLimit',[min(xlon) max(xlon)]);
   h=axesm('MapProjection','eqdcylin','MapLatLimit',[(min(lat)) (max(lat))],'MapLonLimit',[min(lon) max(lon)]);
   pcolorm(lat,lon,varZ); % add additional value to avoid line in pcolor plots
   plotm(mlat,mlon,'color',[0.5 0.5 0.5])
   colormap(cmap) % colormap from 'WhiteYellowOrangeRed'
   caxis([-1 1])
   colorbar
   %set(gca,'fontsize',18)
   title('corr(NASW PDSI, T2m), 10 yrs')
   %set(get(colorbar('location','southoutside'),'xlabel'),'string',sn{j},'fontsize',14)
   %set(get(colorbar('location','southoutside'),'xlabel'),'string','Correlation','fontsize',18)
   %set(cbh,'YTick',[-1:0.2:1])
   %delete( cbh )
   hold off

   %print(['figs/corrmap_naswpdsi_t2m_phyda.png'],'-dpng','-r300');


end










