# Function to plot the (90%) central region using ED
# Inputs: fmat = p x n matrix of curves (each column is a sample),
#         xvals = list of observed points (on x-axis), ie observed wavelengths,
#         ggprint = immediately print plot ('T' or 'F'), where 'F' will save it as a ggplot object
plot_cr = function(fmat, xvals, alpha = 0.9, xlabel = "Wavelength (nm)", ylabel = "% Reflectivity", main = "", ggprint = F) {
  fmat.ed = edepth_set(fmat)
  fmat.cr = central_region(fmat, fmat.ed, alpha)
  fmat.gg = melt(fmat) #Var1 is the wavelength, Var2 is the subj index
  fmat.gg$Var1<-rep(xvals,ncol(fmat))
  cr.gg = melt(cbind(fmat.cr[[1]], fmat.cr[[2]]))
  iq.range = subset(cr.gg,Var2==2)$value-subset(cr.gg,Var2==1)$value
  cr.gg$Var1<-xvals

  med.gg = melt(fmat[,which(fmat.ed == 1)]) #id median curve
  med.gg$Var1 = xvals

  cr.plot = ggplot() +
    geom_point(data = fmat.gg,
               aes(x = Var1,
                   y = value,
                   group = Var2),
               color = "#CCCCCC",
               alpha = 0.3,
               size = 0.05)
  cr.plot<- cr.plot+  geom_point(data = cr.gg,
                                 aes(x = Var1,
                                     y = value,
                                     group = Var2),
                                 color = "#333333",
                                 size = 0.2) +
    geom_point(data = med.gg,
               aes(x = Var1,
                   y = value),
               color = "#CC0000",
               size = 0.2) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(main)

  if (ggprint == T) print(cr.plot)
  return(list(plot=cr.plot))
}


# Function to compare the (90%) central region of 2 datasets
# Inputs: fmat1 = p1 x n1 matrix of curves for dataset 1 (each column is a sample),
#         xvals1 = list of observed points (on x-axis) for dataset 1,
#         fmat2 = p2 x n2 matrix of curves for dataset 1 (each column is a sample),
#         xvals2 = list of observed points (on x-axis) for dataset 2,
#         ggprint = immediately print plot ('T' or 'F'), where 'F' will save it as a ggplot object
#         plot.data = plot the individual sample curves ('T' or 'F')
plot_cr_comp = function(fmat1, fmat2, xvals1, xvals2, alpha = 0.9, xlabel = "Time", ylabel = "% Reflectivity", main = "", ggprint = F,
                        plot.data=T) {
  fmat.ed1 = edepth_set(fmat1)
  fmat.ed2 = edepth_set(fmat2)

  fmat.cr1 = central_region(fmat1, fmat.ed1, alpha)
  fmat.cr2 = central_region(fmat2, fmat.ed2, alpha)

  fmat.gg1 = melt(fmat1) #Var1 is the wavelength, Var2 is the subj index
  fmat.gg1$Var1<-rep(xvals1,ncol(fmat1))
  cr.gg1 = melt(cbind(fmat.cr1[[1]], fmat.cr1[[2]]))
  iq.range1 = subset(cr.gg1,Var2==2)$value-subset(cr.gg1,Var2==1)$value
  cr.gg1$Var1<-xvals1

  med.gg1 = melt(fmat1[,which(fmat.ed1 == 1)])
  med.gg1$Var1 = xvals1

  fmat.gg2 = melt(fmat2) #Var1 is the wavelength, Var2 is the subj index
  fmat.gg2$Var1<-rep(xvals2,ncol(fmat2))
  cr.gg2 = melt(cbind(fmat.cr2[[1]], fmat.cr2[[2]]))
  cr.gg2$Var1<-xvals2

  med.gg2 = melt(fmat2[,which(fmat.ed2 == 1)])
  med.gg2$Var1 = xvals2


  cr.plot<- ggplot()
  if(plot.data==T){
    cr.plot<- cr.plot + geom_point(data = fmat.gg1, #first dataset
                                   aes(x = Var1,
                                       y = value,
                                       group = Var2),
                                   color = "#FFB6C1",
                                   alpha = 0.3,
                                   size = 0.05) +
      geom_point(data = fmat.gg2, #second dataset
                 aes(x = Var1,
                     y = value,
                     group = Var2),
                 color = "#98F5FF",
                 alpha = 0.3,
                 size = 0.05)
  }
  cr.plot<- cr.plot+
    geom_point(data = cr.gg1,
               aes(x = Var1,
                   y = value,
                   group = Var2),
               color = "#FF0000",
               size = 0.5) +
    geom_point(data = med.gg1, #median
               aes(x = Var1,
                   y = value),
               color = "#8B0000",
               size = 0.5) +
    geom_point(data = cr.gg2,
               aes(x = Var1,
                   y = value,
                   group = Var2),
               color = "#0000FF",
               size = 0.5) +
    geom_point(data = med.gg2,
               aes(x = Var1,
                   y = value),
               color = "#00008B",
               size = 0.5) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(main)

  if (ggprint == T) print(cr.plot)
  return(list(plot=cr.plot))
}

# Function to plot a functional box plot using ED
# Inputs: fmat = p x n matrix of curves (each column is a sample),
#         xvals = list of observed points (on x-axis), ie observed wavelengths,
#         ggprint = immediately print plot ('T' or 'F'), where 'F' will save it as a ggplot object
plot_fbplot = function(fmat, xvals, xlabel = "Wavelength (nm)", ylabel = "% Reflectivity", main = "", ggprint = F) {
  fmat.ed = edepth_set(fmat)
  fmat.cr = central_region(fmat, fmat.ed, 0.5)
  fmat.gg = melt(fmat) #Var1 is the wavelength, Var2 is the subj index
  fmat.gg$Var1<-rep(xvals,ncol(fmat))
  cr.gg = melt(cbind(fmat.cr[[1]], fmat.cr[[2]]))
  iq.range = subset(cr.gg,Var2==2)$value-subset(cr.gg,Var2==1)$value
  cr.gg$out.lim<-c(subset(cr.gg,Var2==1)$value-iq.range*1.5,subset(cr.gg,Var2==2)$value+iq.range*1.5) #1 is low limits, #2 is high limit
  cr.gg$Var1<-xvals
  outliers<-out.ind<-NULL
  for(i in unique(fmat.gg$Var2)){
    subj<-subset(fmat.gg,Var2==i)
    check_low<-sum(as.numeric(subj$value<subset(cr.gg,Var2==1)$out.lim))
    check_high<-sum(as.numeric(subj$value>subset(cr.gg,Var2==2)$out.lim))
    if(check_low+check_high>10){
      outliers<-rbind(outliers,subj)
      out.ind<-c(out.ind,i)
    }
  }
  med.gg = melt(fmat[,which(fmat.ed == 1)]) #id median curve
  med.gg$Var1 = xvals

  cr.plot = ggplot() +
    geom_point(data = fmat.gg,
               aes(x = Var1,
                   y = value,
                   group = Var2),
               color = "#CCCCCC",
               alpha = 0.3,
               size = 0.05)
  if(!is.null(outliers)){
    cr.plot<-cr.plot+geom_point(data = outliers,
                                aes(x = Var1,
                                    y = value,
                                    group = Var2),
                                color = "#33FF00",
                                alpha = 0.3,
                                size = 0.1)
  }
  cr.plot<- cr.plot+  geom_point(data = cr.gg,
                                 aes(x = Var1,
                                     y = value,
                                     group = Var2),
                                 color = "#333333",
                                 size = 0.2) +
    geom_point(data = cr.gg,
               aes(x = Var1,
                   y = out.lim,
                   group = Var2),
               color = "#666666",
               size = 0.2) +
    geom_point(data = med.gg,
               aes(x = Var1,
                   y = value),
               color = "#CC0000",
               size = 0.2) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(main)

  if (ggprint == T) print(cr.plot)
  return(list(plot=cr.plot,outlier.ind=out.ind))
}
