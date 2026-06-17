#!/usr/bin/env Rscript
# =====================================================================
# Step 7: transitivity table + prespecified sensitivity analyses.
# PRIMARY = clean network (9 contrasts; confounded head-to-heads Haghighat
# & Gurdal demoted). Re-derives change-score SDs under different pre/post
# correlations and an endpoint scoring, then refits the component model.
# =====================================================================
suppressPackageStartupMessages({library(netmeta); library(metafor); library(dplyr)})
root <- normalizePath(getwd(), mustWork=TRUE)
outdir <- file.path(root,"final_analysis","outputs")
base <- read.csv(file.path(outdir,"cnma_arm_level.csv"), stringsAsFactors=FALSE)

wan_S1 <- function(a,m,b,n){ mean<-(a+2*m+b)/4+(a-2*m+b)/(4*n); sd<-(b-a)/(2*qnorm((n-0.375)/(n+0.25))); c(mean=mean,sd=sd) }
wan_S2 <- function(q1,m,q3,n){ mean<-(q1+m+q3)/3; sd<-(q3-q1)/(2*qnorm((0.75*n-0.125)/(n+0.25))); c(mean=mean,sd=sd) }
sd_change <- function(s1,s2,r) sqrt(s1^2+s2^2-2*r*s1*s2)
CONF <- c("SRMA-0405","SRMA-0905")  # confounded head-to-heads (Haghighat, Gurdal)

# ---- builder: all 11 contrasts given r and scoring -------------------
build_contrasts <- function(r=0.5, scoring=c("change","endpoint")){
  scoring <- match.arg(scoring); rows<-list(); k<-0
  addc <- function(id,t1,t2,n1,m1,s1,n2,m2,s2){ k<<-k+1
    es<-escalc("SMD",m1i=m1,sd1i=s1,n1i=n1,m2i=m2,sd2i=s2,n2i=n2)
    rows[[k]]<<-data.frame(study_id=id,treat1=t1,treat2=t2,TE=as.numeric(es$yi),seTE=sqrt(as.numeric(es$vi))) }
  addc("SRMA-0457","CB+MLD","CB",23,46.1,22.6,21,38.6,16.1)
  se<-(22.9+19)/(2*1.96); sp<-se/sqrt(1/39+1/38); addc("SRMA-0243","CB+MLD","CB",39,1.9,sp,38,0,sp)
  sa<-(65-32)/(2*1.96)*sqrt(20); sc<-(78-43)/(2*1.96)*sqrt(22); addc("SRMA-0901","CB+MLD","CB",20,48,sa,22,60,sc)
  sd_dv<-(1.8+2.1)/(2*1.96)/sqrt(1/64+1/65); addc("SRMA-0144","CB+MLD","CB",64,-0.2,sd_dv,65,0,sd_dv)  # De Vrieze
  addc("SRMA-0469","CB+IPC","CB",12,45.3,18.2,11,26.0,22.1)
  ti<-wan_S2(48.3,54.6,58.2,38); tc<-wan_S2(48.3,49.6,54.9,38); addc("SRMA-0210","CB+IPC","CB",38,ti["mean"],ti["sd"],38,tc["mean"],tc["sd"])
  addc("SRMA-0405","CB+IPC","CB+MLD",56,37.5,14.4,56,43.1,13.7)                    # Haghighat (confounded)
  szp<-(7.93-3.06)/qt(0.975,25)/sqrt(1/14+1/13); addc("SRMA-0407","CB+MLD+IPC","CB+MLD",14,7.93,szp,13,3.06,szp)
  ub<-wan_S1(220,840,3460,15); ue<-wan_S1(60,500,2160,15); cb<-wan_S1(180,630,1820,15); ce<-wan_S1(0,480,1410,15)
  jb_m<-c(579,258); je_m<-c(504,252); jb_s<-c(411,203); je_s<-c(382,193)
  if(scoring=="change"){
    addc("SRMA-0904","CB+IPC","CB",15, ub["mean"]-ue["mean"], sd_change(ub["sd"],ue["sd"],r), 15, cb["mean"]-ce["mean"], sd_change(cb["sd"],ce["sd"],r))
    addc("SRMA-0488","CB+MLD","CB+IPC",12, jb_m[1]-je_m[1], sd_change(jb_m[2],je_m[2],r), 12, jb_s[1]-je_s[1], sd_change(jb_s[2],je_s[2],r))
    addc("SRMA-0905","CB+MLD","CB+IPC",15, 3533-3004, sd_change(739,739,r), 15, 3581-3142, sd_change(783,783,r))  # Gurdal (confounded)
  } else {
    addc("SRMA-0904","CB+IPC","CB",15, ce["mean"], ce["sd"], 15, ue["mean"], ue["sd"])
    addc("SRMA-0488","CB+MLD","CB+IPC",12, je_s[1], je_s[2], 12, je_m[1], je_m[2])
    addc("SRMA-0905","CB+MLD","CB+IPC",15, 3142, 783, 15, 3004, 739)
  }
  bind_rows(rows)
}

fit_cnma <- function(d, label){
  net <- tryCatch(netmeta(TE,seTE,treat1,treat2,study_id,data=d,sm="SMD",common=TRUE,random=TRUE,
                          reference.group="CB",method.tau="REML",details.chkmultiarm=FALSE), error=function(e) NULL)
  if(is.null(net)) return(NULL)
  nc <- tryCatch(netcomb(net, inactive="CB"), error=function(e) NULL)
  if(is.null(nc)){
    cb<-which(net$trts=="CB"); g<-function(t){i<-which(net$trts==t); if(!length(i)) return(c(NA,NA,NA)); c(net$TE.random[i,cb],net$lower.random[i,cb],net$upper.random[i,cb])}
    mld<-g("CB+MLD"); ipc<-g("CB+IPC")
    return(data.frame(analysis=label,k=net$k,MLD=mld[1],MLD_lo=mld[2],MLD_hi=mld[3],IPC=ipc[1],IPC_lo=ipc[2],IPC_hi=ipc[3],
                      combo=NA,combo_lo=NA,combo_hi=NA,tau2=net$tau2,I2=net$I2*100,incoh_p=net$pval.Q.inconsistency))
  }
  mi<-which(nc$comps=="MLD"); ii<-which(nc$comps=="IPC"); cmb<-which(nc$trts=="CB+MLD+IPC"); cb<-which(nc$trts=="CB")
  cb3<-if(length(cmb)==1) c(nc$TE.random[cmb,cb],nc$lower.random[cmb,cb],nc$upper.random[cmb,cb]) else c(NA,NA,NA)
  data.frame(analysis=label,k=net$k,
             MLD=nc$Comp.random[mi],MLD_lo=nc$lower.Comp.random[mi],MLD_hi=nc$upper.Comp.random[mi],
             IPC=nc$Comp.random[ii],IPC_lo=nc$lower.Comp.random[ii],IPC_hi=nc$upper.Comp.random[ii],
             combo=cb3[1],combo_lo=cb3[2],combo_hi=cb3[3],tau2=net$tau2,I2=net$I2*100,incoh_p=net$pval.Q.inconsistency)
}
prim <- function(d) subset(d, !study_id %in% CONF)   # drop confounded head-to-heads

res <- list()
res[[1]] <- fit_cnma(prim(build_contrasts(0.5,"change")), "PRIMARY (clean, r=0.5, change)")
res[[length(res)+1]] <- fit_cnma(build_contrasts(0.5,"change"), "(s1) FULL: add confounded head-to-heads (Haghighat, Gurdal)")
# (a) exclude imputed-variance trials (keep only reported-SD: McNeely, Szuba [Haghighat confounded->out])
imp_ids <- c("SRMA-0243","SRMA-0901","SRMA-0144","SRMA-0904","SRMA-0210","SRMA-0488","SRMA-0407")
res[[length(res)+1]] <- fit_cnma(subset(prim(build_contrasts(0.5,"change")), !study_id %in% imp_ids), "(a) Exclude imputed-variance trials")
# (b) r variation
res[[length(res)+1]] <- fit_cnma(prim(build_contrasts(0.3,"change")), "(b) r = 0.3")
res[[length(res)+1]] <- fit_cnma(prim(build_contrasts(0.7,"change")), "(b) r = 0.7")
# (d) exclude high-RoB (McNeely SRMA-0457 in primary; Gurdal already out)
res[[length(res)+1]] <- fit_cnma(subset(prim(build_contrasts(0.5,"change")), study_id!="SRMA-0457"), "(d) Exclude high-RoB (McNeely)")
# (e) leave-one-out over the primary set
P <- prim(build_contrasts(0.5,"change"))
for(id in unique(P$study_id))
  res[[length(res)+1]] <- fit_cnma(subset(P, study_id!=id), paste0("(e) LOO drop ", id))
# (f) endpoint scoring
res[[length(res)+1]] <- fit_cnma(prim(build_contrasts(0.5,"endpoint")), "(f) Endpoint (not change) scoring")
# (g) exclude De Vrieze (show its stabilising role on MLD)
res[[length(res)+1]] <- fit_cnma(subset(P, study_id!="SRMA-0144"), "(g) Exclude De Vrieze (sham-controlled MLD)")

S <- bind_rows(res); num<-sapply(S,is.numeric); S[num]<-lapply(S[num],function(x) round(x,3))
write.csv(S, file.path(outdir,"cnma_sensitivity.csv"), row.names=FALSE)
cat("=== SENSITIVITY (component SMDs; PRIMARY = clean network) ===\n")
print(S[,c("analysis","k","MLD","MLD_lo","MLD_hi","IPC","IPC_lo","IPC_hi","I2","incoh_p")], row.names=FALSE)

trans <- base %>% transmute(study_id, citation, edge, network_role, phase, outcome_construct,
                            scale_direction, timepoint, n_total=n_index+n_ref, sd_source)
write.csv(trans, file.path(outdir,"transitivity_table.csv"), row.names=FALSE)
cat("\nWrote cnma_sensitivity.csv and transitivity_table.csv\n")
