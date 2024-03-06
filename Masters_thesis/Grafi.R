library(plotly)

b<-plot_ly(type = 'bar') %>%
  add_bars(x=colnames(dgeObj)[1],y = ~dgeObj$samples$lib.size[1]) %>%
  add_bars(x=colnames(dgeObj)[2],y = ~dgeObj$samples$lib.size[2]) %>%
  add_bars(x=colnames(dgeObj)[3],y = ~dgeObj$samples$lib.size[3]) %>%
  add_bars(x=colnames(dgeObj)[4],y = ~dgeObj$samples$lib.size[4]) %>%
  add_bars(x=colnames(dgeObj)[5],y = ~dgeObj$samples$lib.size[5]) %>%
  add_bars(x=colnames(dgeObj)[6],y = ~dgeObj$samples$lib.size[6]) %>%
  layout(yaxis = list(title = 'Velikost knjižnic'), xaxis = list(title='Vzorci'), showlegend=FALSE)
#<-plot_ly(type = 'bar',x=colnames(dgeObj)[1:3],y = ~dgeObj$samples$lib.size[1:3], name='Odzivniki') %>%
#  add_bars(x=colnames(dgeObj)[4:6],y = ~dgeObj$samples$lib.size[4:6], marker=list(color='red'), name='Neodzivniki') %>%
#  layout(yaxis = list(title = 'Velikost knjižnic'), xaxis = list(title='Vzorci'), barmode = 'group')


logcounts1<-data.frame(logcounts)
min_x<-(-1)
max_x<-6
boks1<- plot_ly(type = 'box') %>%
  add_boxplot(y=~logcounts1[,1], name = colnames(logcounts1[1]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts1[,2], name = colnames(logcounts1[2]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts1[,3], name = colnames(logcounts1[3]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts1[,4], name = colnames(logcounts1[4]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts1[,5], name = colnames(logcounts1[5]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts1[,6], name = colnames(logcounts1[6]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  layout(shapes=list(type='line', x0=min_x, x1=max_x, y0=median(logcounts.N), y1=median(logcounts.N)),
         xaxis = list(title = 'Vzorci'), yaxis = list(title = 'Log2(CPM)'), showlegend=FALSE)

logcounts.N1<-data.frame(logcounts.N)
boks2<- plot_ly(type = 'box') %>%
  add_boxplot(y=~logcounts.N1[,1], name = colnames(logcounts.N1[1]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts.N1[,2], name = colnames(logcounts.N1[2]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts.N1[,3], name = colnames(logcounts.N1[3]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts.N1[,4], name = colnames(logcounts.N1[4]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts.N1[,5], name = colnames(logcounts.N1[5]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  add_boxplot(y=~logcounts.N1[,6], name = colnames(logcounts.N1[6]), marker = list(color = 'lightblue'),
              line = list(color = 'lightblue')) %>%
  layout(shapes=list(type='line', x0=min_x, x1=max_x, y0=median(logcounts.N), y1=median(logcounts.N)),
         xaxis = list(title = 'Vzorci'), yaxis = list(title = 'Log2(CPM)'), showlegend=FALSE)

#heatmap
h1<-heatmaply(highly_variable_lcpm, colors = c('blue', 'red'), column_text_angle = 0, xlab='', ylab='', 
              seriate = "mean", col_side_colors = sampleinfo$info,
              plot_method = "plotly")
              
h2<-heatmaply(highly_variable_lcpm.N,colors = c('blue', 'red'), column_text_angle = 0, row_dend_left = TRUE, xlab='', ylab='', 
              seriate = "mean", col_side_colors = sampleinfo$info,
              plot_method = "plotly")


#MDS plot
mds1<-plotMDS(dgeObj,col=col.cell)
m1<-data.frame(mds1$cmdscale.out)
m1$status<-c('odzivnik','odzivnik','odzivnik', 'neodzivnik','neodzivnik', 'neodzivnik')
md1<-plot_ly(x=m1$X1, y=m1$X2, type = 'scatter',mode = 'text', text= row.names(m1), 
             color = sampleinfo$info, colors = c("red", "blue"), name = m1$status)%>%
        layout(xaxis = list(range = c(-2.5,2.5), title='logFC'), yaxis = list(range = c(-2.5,2.5), title='logFC'))

mds2<-plotMDS(dgeObj1,col=col.cell)
m2<-data.frame(mds2$cmdscale.out)
m2$status<-c('odzivnik','odzivnik','odzivnik', 'neodzivnik','neodzivnik', 'neodzivnik')
md2<-plot_ly(x=m2$X1, y=m2$X2, type = 'scatter',mode = 'text', text= row.names(m2), 
             color = sampleinfo$info, colors = c("red", "blue"), name = m2$status)%>%
        layout(xaxis = list(range = c(-2.5,2.5), title='logFC'), yaxis = list(range = c(-2.5,2.5), title='logFC'))


#graf ??  ?
a1<- rowMeans(logcounts1[,-1])
a2<- rowMeans(logcounts1[,-2])
a3<- rowMeans(logcounts1[,-3])
a4<- rowMeans(logcounts1[,-4])
a5<- rowMeans(logcounts1[,-5])
a6<- rowMeans(logcounts1[,-6])


md4<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts1[,1]-a1, type = 'scatter',mode = 'markers', name=sampleinfo$sample[1],color = I("lightblue")) %>%
  layout( xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))
  
md5<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts1[,2]-a2, type = 'scatter',mode = 'markers', name=sampleinfo$sample[2],color = I("lightblue")) %>%
      layout(xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))
            
md6<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts1[,3]-a3, type = 'scatter',mode = 'markers', name=sampleinfo$sample[3],color = I("lightblue")) %>%
      layout(xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))
            
md7<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts1[,4]-a4, type = 'scatter',mode = 'markers', name=sampleinfo$sample[4],color = I("lightblue")) %>%
      layout(xaxis = list(title='logCPM'), yaxis = list(title = 'Log-razmerje'))

md8<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts1[,5]-a5, type = 'scatter',mode = 'markers', name=sampleinfo$sample[5],color = I("lightblue")) %>%
      layout(xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))

md9<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts1[,6]-a6, type = 'scatter',mode = 'markers', name=sampleinfo$sample[6],color = I("lightblue")) %>%
      layout(xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))

md10<- subplot(md4,md5,md6,md7,md8,md9, shareY = TRUE, titleX = TRUE) %>%
          layout(showlegend=FALSE)

for (i in 1:6) {
  mdp1<-plot_ly(type = 'scatter', mode='markers')
  tmp<-plot_ly(x=dgeObj1$AveLogCPM, y=logcounts1[,i]-rowMeans(logcounts1[,-i]), type = 'scatter',mode = 'markers', name=sampleinfo$sample[i]) %>%
    layout(xaxis = list(title='Povprečje log(CPM)'), yaxis = list(title = 'Log-razmerje (this sample vs others)'))
  mdp1<-subplot(mdp1,tmp, shareX = TRUE, shareY = TRUE)
}
for (i in 1:6) {
  print(i)}

for (i in 1:6) {
  mdp2<-plot_ly(type = 'scatter', mode='markers')
  tmp2<-plot_ly(x=dgeObj1$AveLogCPM, y=logcounts.N1[,i]-rowMeans(logcounts.N1[,-i]), type = 'scatter',mode = 'markers', name=sampleinfo$sample[1]) %>%
    layout(title =sampleinfo$sample[i], xaxis = list(title='Povprečje log(CPM)'), yaxis = list(title = 'Log-razmerje'))
  mdp2<-subplot(mdp2,tmp2, shareY = TRUE)
}

a11<- rowMeans(logcounts.N1[,-1])
a12<- rowMeans(logcounts.N1[,-2])
a13<- rowMeans(logcounts.N1[,-3])
a14<- rowMeans(logcounts.N1[,-4])
a15<- rowMeans(logcounts.N1[,-5])
a16<- rowMeans(logcounts.N1[,-6])


md11<-plot_ly(x=dgeObj2$AveLogCPM,y=logcounts.N1[,2]-a12, type = 'scatter',mode = 'markers',color = I("lightblue")) %>%
  layout(xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))

md12<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts.N1[,2]-a12, type = 'scatter',mode = 'markers',color = I("lightblue")) %>%
  layout( xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))

md13<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts.N1[,3]-a13, type = 'scatter',mode = 'markers',color = I("lightblue")) %>%
  layout( xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))

md14<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts.N1[,4]-a14, type = 'scatter',mode = 'markers',color = I("lightblue")) %>%
  layout( xaxis = list(title='logCPM'), yaxis = list(title = 'Log-razmerje'))

md15<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts.N1[,5]-a15, type = 'scatter',mode = 'markers',color = I("lightblue")) %>%
  layout( xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))

md16<-plot_ly(x=dgeObj2$AveLogCPM, y=logcounts.N1[,6]-a16, type = 'scatter',mode = 'markers',color = I("lightblue")) %>%
  layout( xaxis = list(title='logCPM'),yaxis = list(title = 'Log-razmerje'))

md17<- subplot(md11,md12,md13,md14,md15,md16, shareY=TRUE, titleX = TRUE) %>%
  layout(showlegend=FALSE)


png(file="Composition_bias_NORM.png")
par(mfrow=c(1,6))
plotMD(dgeObj1,column = 1)
abline(h=0,col="grey")
plotMD(dgeObj1,column = 2) 
abline(h=0,col="grey")
plotMD(dgeObj1,column = 3)
abline(h=0,col="grey")
plotMD(dgeObj1,column = 4) 
abline(h=0,col="grey")
plotMD(dgeObj1,column = 5)
abline(h=0,col="grey")
plotMD(dgeObj1,column = 6) 
abline(h=0,col="grey")
dev.off()


#edgeR
plotBCV(dgeObj2,xlab="Povprečni log(CPM)", ylab="Biološki koeficient variacije")
b<-cbind(sqrt(dgeObj2$common.dispersion),sqrt(dgeObj2$trended.dispersion),sqrt(dgeObj2$tagwise.dispersion),dgeObj2$AveLogCPM)
b<-data.frame(b)
names(b)<-c('common.dispersion','trended.dispersion', 'tagwise.dispersion', 'AveLogCPM')

x_tick=list(title = 'Povprečje logCPM', autotick=FALSE, ticks = "outside", tick0 = 0, dtick = 5, ticklen = 15)
y_tick=list(title = 'Biološki koeficient variacije', autotick=FALSE, ticks = "outside", tick0 = 0, dtick = 0.2, ticklen = 1.4)
bcv<-plot_ly(type = 'scatter', mode ='markers')%>%
  add_markers(x=b$AveLogCPM, y=b$tagwise.dispersion, marker = list(color = 'black',size = 3),name = 'Tagwise') %>%
  add_lines(x=b$AveLogCPM, y=b$common.dispersion, name = 'Common') %>%
  add_lines(x=b$AveLogCPM, y=~b$trended.dispersion, name = 'Trend') %>%
  layout(xaxis = x_tick, yaxis = y_tick)


#p<-ggplot(data=Rezultati, aes(x=limma_logFC, y=edgeR_logFC, group=1)) + geom_point(colour="lightblue") + geom_smooth(method = "lm", col="black") 

p<-plot_ly(x=Rezultati$limma_logFC, y=Rezultati$edgeR_logFC, type='scatter', mode='markers', color = I("lightblue"))%>%
            layout(yaxis = list(title = 'edgeR'), xaxis = list(title = 'limma-VOOM'), showlegend = FALSE) %>%
              add_trace(x=Rezultati$limma_logFC[7322], y=Rezultati$edgeR_logFC[7322], marker = list(color = I("red"))) %>%
              add_markers(x=Rezultati$limma_logFC[1366], y=Rezultati$edgeR_logFC[1366], marker = list(color = I("red"))) %>%
              add_markers(x=Rezultati$limma_logFC[13211], y=Rezultati$edgeR_logFC[13211], marker = list(color = I("red"))) %>%
              add_markers(x=Rezultati$limma_logFC[13193], y=Rezultati$edgeR_logFC[13193], marker = list(color = I("red")))%>%
              add_annotations(x = Rezultati$limma_logFC[c(7322,1366,13211,13193)], y = Rezultati$edgeR_logFC[c(7322,1366,13211,13193)], text = Rezultati$SYMBOL[c(7322,1366,13211,13193)],
                              showarrow = TRUE,  arrowhead = 0, arrowsize = 0, ax = 0,ay =-10)#%>%
              #add_lines(x =Rezultati$limma_logFC, y = fitted(lm( edgeR_logFC ~ limma_logFC, data = Rezultati)), color = I("black"))

#pl<-ggplotly(p)
pl<-style(p,text = paste('</br> Gen: ', Rezultati$SYMBOL,
                          '</br> edgeR: ', Rezultati$edgeR_logFC,
                          '</br> VOOM: ', Rezultati$limma_logFC),
          hoverinfo = 'text')

p1<-plot_ly(x=Rezultati$limma_PValue, y=Rezultati$edgeR_PValue, type='scatter', mode='markers', color = I("lightblue"))%>%
        layout(yaxis = list(title = 'edgeR'), xaxis = list(title = 'limma-VOOM'), showlegend = FALSE) %>%
          add_trace(x=Rezultati$limma_PValue[7322], y=Rezultati$edgeR_PValue[7322], marker = list(color = I("red"))) %>%
          add_markers(x=Rezultati$limma_PValue[1366], y=Rezultati$edgeR_PValue[1366], marker = list(color = I("red"))) %>%
          add_markers(x=Rezultati$limma_PValue[13211], y=Rezultati$edgeR_PValue[13211], marker = list(color = I("red"))) %>%
          add_markers(x=Rezultati$limma_PValue[13193], y=Rezultati$edgeR_PValue[13193], marker = list(color = I("red")))%>%
          add_annotations(x = Rezultati$limma_PValue[c(7322,1366,13211,13193)], y = Rezultati$edgeR_PValue[c(7322,1366,13211,13193)], text = Rezultati$SYMBOL[c(7322,1366,13211,13193)],
                          showarrow = TRUE,  arrowhead = 0, arrowsize = 0, ax = 0,ay =-10)

#p1<-ggplot(data=Rezultati, aes(x=limma_PValue, y=edgeR_PValue, group=1)) + 
  geom_point(colour="lightblue", size=1) + 
  labs(title="p-vrednosti",x="VOOM", y="edgeR")+
  geom_point(data=Rezultati, aes (x=limma_PValue[7322], y=edgeR_PValue[7322], colour = "red"),size=1)+
  geom_point(data=Rezultati, aes (x=limma_PValue[1366], y=edgeR_PValue[1366], colour = "red"),size=1)+
  geom_point(data=Rezultati, aes (x=limma_PValue[13211], y=edgeR_PValue[13211], colour = "red"),size=1)+
  geom_point(data=Rezultati, aes (x=limma_PValue[13193], y=edgeR_PValue[13193], colour = "red"),size=1)+
  theme(legend.position="none")
#pl1<-ggplotly(p1)
  
pl1<-style(p1,text = paste('</br> Gen: ', Rezultati$SYMBOL,
                            '</br> edgeR: ', Rezultati$edgeR_PValue,
                            '</br> VOOM: ', Rezultati$limma_PValue),
           hoverinfo = 'text')


#VOOM
y<-voom(dgeObj1,design,plot=TRUE, save.plot = TRUE) 
summa.fit <- decideTests(fit)
summary(summa.fit)


v<-plot_ly(type='scatter', mode='markers')%>%
  add_trace(x=y$voom.xy$x, y=y$voom.xy$y, marker = list(color = I("black"), size = 4))%>%
  add_lines(x=y$voom.line$x, y=,y$voom.line$y, color = I("red")) %>%
  layout(yaxis = list( title = 'sqrt(standardni odklon)'), xaxis = list(range = c(0,20),title = 'logCPM'), showlegend = FALSE)


mds3<-plotMDS(y)
m3<-data.frame(mds3$cmdscale.out)
m3$status<-c('odzivnik','odzivnik','odzivnik', 'neodzivnik','neodzivnik', 'neodzivnik')
md3<-plot_ly(x=m3$X1, y=m3$X2, type = 'scatter',mode = 'text', text= row.names(m3), 
             color = sampleinfo$info, colors = c("red", "blue"), name=m3$status)%>%
  layout(xaxis = list(range = c(-2.5, 2.5), title='logFC'), yaxis = list(range = c(-2.5, 2.5), title = 'logFC'))





