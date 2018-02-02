setwd('~/applications/Rnall2.0')

d = read.csv('test2.out', header=F, sep='\t', skip=1)
colnames(d) = c('window_size','center_position','delta_g')

d[d$delta_g > 0, 'delta_g'] = 0

d1 = d[d$window_size==20,]
d2 = d[d$window_size==32,]
d3 = d[d$window_size==50,]
d0 = d1
d0$delta_g = 0

plot(d3$center_position, d3$delta_g, type='l', lwd=1, col='#000000'
     , main='Spooky energy landscape\n(local structures)', xlab='POSITION', ylab='d.G'
     , xaxs='i', yaxs='i')
points(d2$center_position, d2$delta_g, type='l', lwd=1, col='#000000')
points(d1$center_position, d1$delta_g, type='l', lwd=1, col='#FF0000')

polygon(c(d3$center_position,rev(d3$center_position)), c(rep(0, nrow(d3)), rev(d3$delta_g)), col='#FF0000', border=NA)
polygon(c(d2$center_position,rev(d2$center_position)), c(rep(0, nrow(d2)), rev(d2$delta_g)), col='#770000')
polygon(c(d1$center_position,rev(d1$center_position)), c(rep(0, nrow(d1)), rev(d1$delta_g)), col='#220000')

legend('bottom',c('20-bases','32-bases','50-bases')
       , fill=c('#220000','#770000','#FF0000'), title='Sequence length')