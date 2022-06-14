#!usr/bin/R

# R to IDL script

suppressMessages(library(nlme))
suppressMessages(library(mgcv))
suppressMessages(library(segmented))
suppressMessages(library(pracma))
suppressMessages(library(spatstat))

max_interpol <- function(y,x,xlarge,inds){
  outy = vector(mode='numeric', length = length(xlarge))
  outy[outy == 0] = NA
  #get rid of na values in y
  goody = y[!is.na(y)]
  goodx = x[!is.na(y)]
  mpinds = inds[!is.na(y)]
  
  #mpinds = which(xlarge %in% goodx) # this didn't work sometimes multiple vals of the same dist if SIC is 0
  outy[mpinds] = goody
  
  if (length(outy[mpinds])!=length(goody)) browser()
  
  #novalinds = which(!(xlarge %in% goodx))
  novalinds = which(is.na(outy)) # find remaining vals to fill
  
  for (i in 1:length(novalinds)){
    #find closest value to left and right
    thisind = novalinds[i]
    leftind = which(mpinds < thisind)
    rightind = which(mpinds > thisind)
    
    if (length(leftind) > 0) mpleftind = max(leftind) else mpleftind = 1
    if (length(rightind) > 0) mprightind = min(rightind) else mprightind = length(mpinds)
    
    outy[thisind] =  max(c(goody[mpleftind],goody[mprightind]))
  }
  return(outy)

}

args = c('/jill_hdd/Batch_process/2019_12/mizfind/csv',
        '/jill_hdd/Batch_process/2019_12/mizfind/rawheights',
         '/jill_hdd/Batch_process/2019_12/mizfind/plots',
         '/jill_hdd/Batch_process/2019_12/mizfind/Rout')

#args = c('/jill_hdd/Batch_process/2019_02/mizfind/csv',
#         '/jill_hdd/Batch_process/2019_02/mizfind/rawheights',
#         '/jill_hdd/Batch_process/2019_02/mizfind/plots',
#         '/jill_hdd/Batch_process/2019_02/mizfind/Rout')

#args <- commandArgs(trailingOnly=TRUE)
# want the csv 

csvstring = args[1]
rawfilestring = args[2]
plotout = args[3]
out = args[4]
setwd(csvstring)

filenames = list.files(path = csvstring, pattern = '*.csv')
rawheightfiles = list.files(path=rawfilestring, pattern = '*.csv')

# set up df for saving good filenames as csv later
goodfiles = data.frame(matrix(ncol=14*2+1,nrow=length(filenames)))

d = read.csv(filenames[1], header=T)
col=(c('LineID', paste(names(d[3:16]),'_DC',sep=''), paste(names(d[3:16]),'_D',sep='')))
colnames(goodfiles) = col
goodfilecount = 1

dmax = 1000



#filenames = filenames[grep('*20190201093627_05340201_end-*', filenames)]

for (f in 187:length(filenames)){
  
  print(paste(f,'out of',length(filenames)))
  
  d <- read.csv(filenames[f], header=T)
  lineID = strsplit(filenames[f],'-')[[1]][1]
  
  # open pdf for saving plot
  
  plotpath= paste(plotout,'/',lineID,'-mizcheck.pdf', sep='')
  pdf(file=plotpath)
  
  # get matching raw height values and plot them
  
  rawheightmatch = grep(lineID, rawheightfiles)
  rawheightd = read.csv(paste(rawfilestring,'/',rawheightfiles[rawheightmatch],sep=''), header=T)

  par(mfrow = c(3, 1))
  par(mgp=c(2,1,0))
  par(mar = c(3, 3, 2, 1))
  par(xpd=FALSE)
  
  b1dist = rawheightd$b1_dist/1000
  b2dist = rawheightd$b2_dist/1000
  b3dist = rawheightd$b3_dist/1000
  
  ylimmin = max(c(min(rawheightd[,c(2,4,6)],na.rm=T),-5))
  ylimmax = min(c(max(rawheightd[,c(2,4,6)],na.rm=T),5))
  
  plot(rawheightd$b1_height~b1dist, type="l", lty=1,ylim=range(ylimmin,ylimmax),xlab='Along-track distance (km)',ylab = 'Height (m)',main = 'Beam 1')
  plot(rawheightd$b2_height~b2dist,type='l',lty=1,ylim=range(ylimmin,ylimmax),xlab='Along-track distance (km)',ylab = 'Height (m)',main = "Beam 2")
  plot(rawheightd$b3_height~b3dist,type='l',lty=1,ylim=range(ylimmin,ylimmax),xlab = 'Along-track distance (km)', ylab = 'Height(m)',main = 'Beam 3')
  
  
   # reference is [row,column]
  
  # set up df for storing fit results
  #fit params (slope/int) and errors
  #break point (bp) index
  #MIZ estimate - where left fitted break (attenuation part) intersects quad error
  #RMSE of fit - predicted line vs actual values within miz dist
  #error of fit - sum of relative errors in slope and intercept * by distance, note: only within segmented bounds not whole miz dist
  
  fitresd = data.frame(matrix(nrow=14, ncol=49))
  colnames(fitresd) = c('DC_slope_lin','DC_slope_lin_err',
                        'DC_int_lin','DC_int_lin_err',
                        'DC_bp_lin', 'DC_miz_lin', 
                        'DC_RMSE_lin','DC_miz_lin_err', # 1-8
                        
                        'DC_slope_exp','DC_slope_exp_err',
                        'DC_exp_lin','DC_int_exp_err',
                        'DC_bp_exp', 'DC_miz_exp', 
                        'DC_RMSE_exp','DC_miz_exp_err', #9-16
                        
                        'D_slope_lin','D_slope_lin_err',
                        'D_int_lin','D_int_lin_err',
                        'D_bp_lin', 'D_miz_lin', 
                        'D_RMSE_lin','D_miz_lin_err', #17:24
                        
                        'D_slope_exp','D_slope_exp_err',
                        'D_exp_lin','D_int_exp_err',
                        'D_bp_exp', 'D_miz_exp', 
                        'D_RMSE_exp','D_miz_exp_err', # 25-32
                        
                        # add in breakpoint info in case we go with this method
                        'DC_segslop_lin','DC_segslop_lin_err', 
                        'DC_segint_lin', 'DC_bp_lin_err', #33-36
                        
                        'DC_segslop_exp','DC_segslop_exp_err',
                        'DC_segint_exp', 'DC_bp_exp_err', #37-40
                        
                        'D_segslop_lin','D_segslop_lin_err',
                        'D_segint_lin','D_bp_lin_err',#41-44
                        
                        'D_segslop_exp','D_segslop_exp_err',
                        'D_segint_exp','D_bp_exp_err', #45-48
                        
                        'SIC_dist'
                        ) 
  rownames(fitresd) = names(d[3:16])     
  
  mpinds = d$good_mpinds+1 # need the plus one because its in R
  
  for (hsm in 1:14){
    
    # set up new plot page
    #par(mfcol = c(3, 2))
    
    
    layout_matrix = matrix(c(1,2,3,1,4,5),ncol=2,nrow=3)
    layout(mat = layout_matrix,
           heights = c(1, 1,1), # Heights of the two rows
           widths = c(2, 2))
    
    #par(mfcol = c(3, 2))
    par(mgp=c(2,1,0))
    par(mar = c(3, 3, 2, 1))
    par(xpd=FALSE)
    
    mhsind = hsm+2
    errind = mhsind+14
    lab = names(d)[mhsind]
    
    # do gam on dist only so we are fitting to the same data
    thisdf = d[d$Dist<dmax,c(1,2,mhsind,errind)]
    names(thisdf) = c('x','xcorr','y','err')
    weights = 1/thisdf$err
    thisdf$normweights = weights/mean(weights, na.rm=T)
    
    if (length(na.omit(thisdf$y))>10) {
      gamlin <- gam(y~s(x),data=thisdf, weights=thisdf$normweights, na.action=na.omit)
      maxx <- max(thisdf$x, na.rm=T)
      minx <- min(thisdf$x, na.rm=T)
      d.pr <- data.frame(x=d$all_SIC_dist)
      d.pr$Pr <- predict(gamlin, newdata=d.pr, type = 'response')
      
      minpoints = findpeaks(as.vector(-(d.pr$Pr)), nups=1, ndowns=1, zero = "-")
      firstpeakind = minpoints[1,2]
      
      # plot the GAM
      plot(y~x, data=thisdf, xlab = 'Distance from ice edge (km)', ylab = 'Mean Hs (m)')
      lines(Pr~x, data=d.pr, col='Dodger blue')
      mtext(paste(lineID, 'method',lab, sep = ' '), side=3, cex = 0.7)
    } else firstpeakind = NA

    
    #GAM criteria are local min 
    if (length(firstpeakind) > 0) {
      peakdist = d.pr$x[firstpeakind] 
      snipind = which(thisdf$x==max(thisdf$x[thisdf$x<=2*peakdist],na.rm=T))
      peakdistind = which(thisdf$x==max(thisdf$x[thisdf$x<=peakdist],na.rm=T))
      if (is.na(firstpeakind)==T) peakdist = -1
    }else {
      peakdist = -1
      snipind = NA
    }
    
    if (peakdist > 0) {
      
      abline(v=peakdist, col='firebrick')
      
      for (distm in 1:2){
        
        if (distm == 1) {
          #thisdf = d[d$Dist<dmax,c(2,mhsind,errind)]
          alldist = d$all_SIC_corrdist
          xlab = 'Corrected Distance from ice edge (km)'
          distind =2
          peakdist = thisdf$xcorr[peakdistind]
          #peakdist = d$all_SIC_corrdist[firstpeakind]
        }
        if (distm == 2) {
          #thisdf = d[d$Dist<dmax,c(1,mhsind,errind)]
          alldist = d$all_SIC_dist
          xlab = 'Distance from ice edge (km)'
          distind =1
          peakdist = thisdf$x[peakdistind]
          #peakdist = d$all_SIC_dist[firstpeakind]
        }

        # broken regression
        estseg = peakdist
        #snipdf = thisdf[thisdf[distind]<2*peakdist,]
        snipdf = thisdf[0:snipind,]
        x= snipdf[,distind]
        y= snipdf[,3]
        err = snipdf$normweights
        

        # see if there are 5 obs before and after
        numbefore = length(na.omit(x[x<peakdist]))
        numafter = length(na.omit(x[x>peakdist]))
        
        # 7 is the lower limit of points needed on either side of break for the psi thing to work
        if ((numbefore >10)&(numafter>10)){
          
          fit0 <- lm(y ~ x, weights=err)
          seglin <- segmented(fit0,psi=estseg,npsi=1,seg.Z=~x) # weight by err and tell it that x is the response
          seglinsum = summary(seglin, var.diff = T, p.df=0) # use diff error variance estimates per section and divide RSS by mean number of params per segment (2)
          
          
          # extract params for lines
          linbreak = seglinsum$psi[[2]]
          linslop = seglinsum$Ttable[[2]]
          linint = seglinsum$Ttable[[1]]
          
          
          # calclate stuff for mizfinding
          if (length(linbreak) >0) {
            secslop = coef(seglin)[3]+linslop
            secint = linint+linslop*linbreak-secslop*linbreak # solve for height of second line
            interpolerr = max_interpol(d[,errind], d[,distind], alldist,mpinds)
            
            if (length(interpolerr) != length(alldist)) browser()
            
            # calculate maximum fitted value 
            #mindist = min(x, na.rm=T)
            #maxlinval = linslop*mindist+linint
            
            # criteria for selection: negative slope and more than 50% data in dist est
            
            lingoodflag = 1
            if (length(linbreak)<1) lingoodflag=0
            if (length(linbreak)>1) if (linbreak <=0) lingoodflag = 0
            if (linslop > 0) lingoodflag = 0 
            
            #npoints = length(na.omit(x[x<linbreak]))
            #if ((linbreak > 0) & (npoints < linbreak/2)) lingoodflag = 0
            #if ((linbreak >0)&(maxlinval > max(thisdf$err[thisdf$x<linbreak], na.rm=T))) lingoodflag = 1 else lingoodflag =0
            #if ((linbreak >0)&(maxlinval<0.25)) lingoodflag =0 
            
            # calculate lines
            linpred = linslop*alldist+linint
            secpred = secslop*alldist+secint
            
            # check mizend - get rid if mizend isn't within the data range
            mizend = min(which(linpred<interpolerr))
            if (is.infinite(mizend)==T) lingoodflag = 0
            #if ((is.infinite(mizend)==F) & (alldist[mizend] > 2*peakdist)) lingoodflag = 0
            #if (mizend > length(na.omit(d$all_SIC_dist))) lingoodflag = 0 # interpolerr is only for this range so this criteria is obsolete
            if (is.infinite(mizend)==F) if (alldist[mizend] < min(x,na.rm=T)) lingoodflag=0
            if (lingoodflag ==1) goodlab = 'Good' else goodlab = 'Bad'
            
            if ((is.infinite(mizend)==F) & (alldist[mizend] > max(x, na.rm=T))) {
              plot(thisdf$y~thisdf[,distind],xlab=xlab,ylab='mean Hs (m)',ylim=range(0,max(thisdf$y,na.rm=T)))
            } else plot(y~x, xlab=xlab, ylab = 'Mean Hs (m)',ylim=range(0,max(y,na.rm=T)))
            
            lines(linpred~alldist)
            lines(secpred~alldist)
            abline(v=linbreak,  lty=3)
            lines(interpolerr~alldist, col = 'dark gray')
            mtext(goodlab, side = 3)
            
            
            if (lingoodflag == 1){
              
              # calculate error params
              mizdist = d$all_SIC_dist[mizend] 
              abline(v=alldist[mizend], col='red')
              
              # calculate RMSE of fit 
              fitdist = alldist[mizend]
              fity = y[x<fitdist]
              fitx = x[x<fitdist]
              mizlinpred = linslop*fitx+linint
              RMSElin = sqrt(mean((mizlinpred-fity)^2, na.rm=T))
              
              # calculate error in MIZ est
              mizerr = mizdist*(seglinsum$Ttable[[5]]/abs(linslop)+seglinsum$Ttable[[4]]/linint)
              
              all_lin_res = list(linslop, seglinsum$Ttable[[5]], # lin slope and std error
                                 linint, seglinsum$Ttable[[4]],  # lin int and std error
                                 linbreak,  mizdist, RMSElin, mizerr)
              extravals = list(secslop,seglinsum$Ttable[3,2],secint,seglinsum$psi[[3]])
              if (distm==1) {
                fitresd[hsm,1:8] = all_lin_res
                fitresd[hsm,33:36] = extravals
              }
              if (distm==2) {
                fitresd[hsm,17:24] = all_lin_res
                fitresd[hsm,41:44] = extravals
              }
              
            } # lingoodflag
          }# breakpoint defined
          
          # same for exponential fit 
          fit0 <- lm(log(y) ~ x, weights=err)
          segexp <- segmented(fit0,psi=estseg,npsi=1,seg.Z=~x) # weight by err and tell it that x is the response
          segexpsum = summary(segexp, var.diff = T, p.df=0)
          
          expbreak = segexpsum$psi[[2]]
          expslop = segexpsum$Ttable[[2]]
          expint = segexpsum$Ttable[[1]]
          
          # calclate stuff for mizfinding
          secslop = coef(segexp)[3]+expslop
          secint = expint+expslop*expbreak-secslop*expbreak # solve for height of second line
          interpolerrexp = log(interpolerr)

          if (length(expbreak) >0) {
            # criteria for selection: negative slope and more than 50% data in dist est
            expgoodflag = 1
            if (length(expbreak)<1) expgoodflag=0
            if (length(expbreak)>1) if (expbreak <=0) expgoodflag = 0
            if (expslop > 0) expgoodflag = 0 
            #npoints = length(na.omit(x[x<expbreak]))
            #if ((expbreak > 0) & (npoints < expbreak/2)) expgoodflag = 0
            #if ((expbreak >0)&(maxlinval > max(thisdf$err[thisdf$x<linbreak], na.rm=T))) lingoodflag = 1 else lingoodflag =0
            #if ((expbreak >0)&(maxexpval<0.25)) expgoodflag =0 
            
            exppred = expslop*alldist+expint
            secpred = secslop*alldist+secint
            
            mizend = min(which(exppred<interpolerrexp))
            if (is.infinite(mizend)==T) expgoodflag = 0
            if (is.infinite(mizend)==F) if (alldist[mizend] < min(x,na.rm=T)) expgoodflag=0
            #if ((is.infinite(mizend)==F) & (alldist[mizend] > max(x,na.rm=T))) expgoodflag = 0
            #if ((is.infinite(mizend)==F) & (alldist[mizend] > 2*peakdist)) expgoodflag = 0
            #if (mizend > length(na.omit(d$all_SIC_dist))) expgoodflag = 0
            if (length(interpolerrexp) != length(exppred)) browser()
            
            if (expgoodflag ==1) goodlab = 'Good' else goodlab = 'Bad'
 
          
            if ((is.infinite(mizend)==F) & (alldist[mizend] > max(x, na.rm=T))) {
              plot(log(thisdf$y)~thisdf[,distind],xlab=xlab,ylab='Log(mean Hs (m))',ylim=range(min(log(interpolerr)),max(log(thisdf$y),na.rm=T)))
            } else plot(log(y)~x, xlab=xlab, ylab = 'Log(Mean Hs (m))',ylim=range(min(log(interpolerr)),max(y,na.rm=T)))
            
            # make plot and label with goodflag
            lines(exppred~alldist)
            lines(secpred~alldist)
            abline(v=expbreak,  lty=3)
            lines(interpolerrexp~alldist, col = 'dark gray')
            mtext(goodlab, side = 3)
            
             
            if (expgoodflag == 1){
              mizdist = d$all_SIC_dist[mizend] 
              abline(v=alldist[mizend], col='red')
              
              # calculate RMSE of fit 
              fitdist = alldist[mizend]
              fityexp = y[x<fitdist]
              fitxexp = x[x<fitdist]
              mizexppred = exp(expslop*fitxexp+expint)
              RMSEexp = sqrt(mean((mizexppred-fityexp)^2, na.rm=T))
              
              # calculate error in MIZ est
              mizerr = mizdist*(segexpsum$Ttable[[5]]/abs(expslop)+segexpsum$Ttable[[4]]/expint)
              
              all_exp_res = list(expslop, segexpsum$Ttable[[5]], # lin slope and std error
                                 expint, segexpsum$Ttable[[4]],  # lin int and std error
                                 expbreak,  mizdist, RMSEexp, mizerr)
              extravals_exp =  list(secslop,segexpsum$Ttable[3,2],secint,segexpsum$psi[[3]])
              if (distm==1) {
                fitresd[hsm,9:16] = all_exp_res
                fitresd[hsm,37:40] = extravals_exp
              }
              if (distm==2) {
                fitresd[hsm,25:32] = all_exp_res
                fitresd[hsm,45:48] = extravals_exp
              }

              
            } #expgoodflag
          } #expbreak
        } # estseg
        
    }#distm

    }#peakdist
      

  } #hsm
  
  # get sic dist
  SICmizdistind = which(d$all_SIC_dist == min(d$all_SIC_dist[d$all_SIC_vals > 80],na.rm=T))
  SICmizdist = d$all_SIC_dist[SICmizdistind]
  fitresd[1,49]=SICmizdist
  
  dev.off()

  # if there are some results save the fitresults as csv per line
  natest = fitresd[1:11,] # don't use the last three firfs they're shit!
  numna = sum(is.na(natest))
  numtot = ncol(natest) * nrow(natest)
  if (numna < numtot) {
    outpath = paste(out,'/', lineID,'.csv', sep='')
    write.csv(fitresd, outpath)
    goodfiles[goodfilecount,1] = lineID
    goodfilecount = goodfilecount+1
  } else file.rename(plotpath, paste(plotout,'/bad/',paste(lineID,'-mizcheck.pdf',sep=''),sep=''))
  
}






